// File: BasisSet/Imp/DeltaFit_IBS.C  The delta basis's collocation tables.
//
// Phi_gi = chi_i(r_g): the values of an ORBITAL block's functions at MY points.  It lives here, and not in
// the XC quadrature that consumes it, because it is an integral over this basis's own functions -- the
// R1.0 discipline ("code that wants an integral asks the basis for it") applied to the one table the
// DENSITY needs in order to contract its private D against this quadrature.  Moving it here is what let
// DM_RhoAtPoints stop taking a raw point list, i.e. what closed the last coordinate escape.
//
// Cached per BLOCK and geometry-fixed: built once per run per block, serving every SCF iteration.  Keyed
// by the block's spatial Irrep (never a pointer), with the real and complex tables in separate maps --
// their key sets are disjoint, because a real TRIM block is its own irrep (doc/RealComplexPlan.md 3c-3).
module;
#include <algorithm>
#include <cassert>
#include <cmath>       // std::abs (the Phi-sparsity magnitude screen)
#include <complex>     // std::abs on a complex table entry
#include <map>         // the per-block table caches
#include <exception>
#include <iostream>
#include <vector>
module qchem.BasisSet.DeltaFit_IBS;
import qchem.Blaze;
import qchem.Reporting;   // report::Timed -- the setup timing bucket
import qchem.Parallel;    // WorkerThreads (GPW_OMP_THREADS): the table build is per-point-block parallel

namespace qchem::BasisSet
{

// GPW_PHI_SPARSITY=1: how sparse is Phi REALLY?  The Phi-SPARSITY item (doc/OpenWork.md Step 3) asserts
// that Phi is stored dense (npts x n) while the true object is sparse, so batching the mesh with a
// per-batch SIGNIFICANT-FUNCTION list makes every Phi-shaped cost -- the rho GEMM and the H_xc GEMM --
// O(npts n_sig^2).  That is a claim about a number nobody has measured, and the win is quadratic in it,
// so measure before building the machinery (the standing lesson of this campaign).
//
// Two granularities, because the mesh already carries BOTH batchings and neither had to be built:
//   SITE   -- the per-atom blocks the Becke build records (Mesh::NSites/SiteBegin/SiteEnd);
//   BATCH  -- contiguous runs of points, which INSIDE a site are radial-shell groups, because a site's
//             points are laid out radial-outer x angular-inner by ProductMesh.
// The reported ratio is the CEILING on the GEMM win: sum_b npts_b n_sig_b^2 against npts n^2.
template <class U> static void ReportPhiSparsity(const mat_t<U>& P, const qcMesh::Mesh& mesh)
{
    static const bool on=std::getenv("GPW_PHI_SPARSITY")!=nullptr;
    if (!on) return;
    const size_t npts=P.rows(), n=P.columns();
    if (npts==0 || n==0) return;
    const double eps=1e-10;                       // the house pointwise magnitude floor (PGData::Reaches)

    auto sigIn=[&](size_t lo, size_t hi)          // functions with ANY significant value on [lo,hi)
    {
        size_t c=0;
        for (size_t i=0;i<n;i++)
            for (size_t q=lo;q<hi;q++)
                if (std::abs(P(q,i))>eps) { c++; break; }
        return c;
    };
    size_t nnz=0;
    for (size_t q=0;q<npts;q++) for (size_t i=0;i<n;i++) if (std::abs(P(q,i))>eps) nnz++;

    const double dense=double(npts)*double(n)*double(n);
    double wSite=0.0; size_t sigMax=0; double sigMean=0.0;
    const size_t nsite=mesh.NSites();
    for (size_t a=0;a<nsite;a++)
    {
        const size_t lo=mesh.SiteBegin(a), hi=mesh.SiteEnd(a), sg=sigIn(lo,hi);
        wSite+=double(hi-lo)*double(sg)*double(sg);
        sigMax=std::max(sigMax,sg); sigMean+=double(sg)/double(std::max<size_t>(nsite,1));
    }
    std::cout<<"[Phi sparsity] npts="<<npts<<" n="<<n
             <<"  nonzero="<<100.0*double(nnz)/(double(npts)*double(n))<<"%"<<std::endl;
    if (nsite>0)
        std::cout<<"[Phi sparsity] SITE batching ("<<nsite<<" blocks): n_sig mean="<<sigMean
                 <<" max="<<sigMax<<" of "<<n<<"  GEMM ceiling="<<dense/std::max(wSite,1.0)<<"x"<<std::endl;
    else
        std::cout<<"[Phi sparsity] SITE batching UNAVAILABLE: this mesh carries no site blocks"<<std::endl;
    // Batch-size sweep: the ceiling rises as batches shrink (more locality) but the per-batch bookkeeping
    // and the GEMM's own efficiency fall, so the useful answer is the TREND, not one number.
    for (size_t bs : {64ul, 256ul, 1024ul, 4096ul})
    {
        double w=0.0; size_t sgMax=0;
        for (size_t lo=0; lo<npts; lo+=bs)
        {
            const size_t hi=std::min(npts,lo+bs), sg=sigIn(lo,hi);
            w+=double(hi-lo)*double(sg)*double(sg); sgMax=std::max(sgMax,sg);
        }
        std::cout<<"[Phi sparsity] BATCH "<<bs<<" pts: n_sig max="<<sgMax<<" of "<<n
                 <<"  GEMM ceiling="<<dense/std::max(w,1.0)<<"x"<<std::endl;
    }
}


template <class U> static mat_t<U> MakePhiAt(const Orbital_1E_IBS<U>* bs, const rvec3vec_t& R, const char* bucket)
{
    qchem::report::Timed timed(bucket);
    // THE dominant SETUP bucket on an atom-centred XC mesh (MnO Γ, 48k points x 122 Bloch functions:
    // 56 s serial, and it DOUBLES on an imposed run's invariant mesh).  Each point is an independent
    // Bloch image sum writing its OWN row, so the loop parallelises with no reduction and no ordering
    // question -- the table is bit-identical at any thread count.  Opt-in (GPW_OMP_THREADS) like every
    // other parallel region here; see qchem.Parallel for why serial is the default.
    //
    // FILL LOCALITY: A CLOSED QUESTION -- DO NOT "FIX" THIS (measured 2026-08-20).  `mat_t` is blaze
    // COLUMN-major and the sweep produces a point's whole ROW, so P(g,0..n-1) strides by npts per element,
    // which LOOKS like the classic cache disaster.  It is not, and the reason is worth keeping: element
    // (g,i) sits at i*npts+g, so CONSECUTIVE POINTS write ADJACENT elements of the SAME cache line.  This
    // is 122 interleaved SEQUENTIAL streams, not a scatter -- one miss per 8 points per function, which
    // the prefetchers handle.  Timed in isolation at this exact shape (99370 x 122): 13 ms column-major
    // vs 8 ms row-major, i.e. ~1.1 ns/element, already the memory-BANDWIDTH floor for writing a 97 MB
    // table at all -- and the row-major route then owes a 21 ms transpose, so it is a NET LOSS as well as
    // a transient second copy of the table.  The whole fill is 0.04% of this bucket; the cost that WAS
    // here lived in the per-image spherical transform (see PG_Spherical's LatticeView) and in the missing
    // pointwise magnitude screen (PGData::Reaches), both now fixed.
    // BATCHED, IN CHUNKS.  The basis is asked for a whole BLOCK of points (VectorFunction's point-set
    // op), not one point at a time, because a periodic basis can then do work once per point that the
    // pointwise call forces it to do once per lattice IMAGE -- see LatticeSum1E::BlochPointValues.  The
    // default op just loops pointwise, so a basis that does not override is unaffected.
    // Chunked rather than one call for the whole mesh so THIS loop keeps the threading: qcMath is a leaf
    // and cannot host OpenMP, and molecular bases reach here too (CreateXCQuadrature is a generic basis
    // face), so handing the whole set to a serial default would deserialise them.  Each chunk writes its
    // own disjoint rows, so the table stays bit-identical at any thread count.
    const size_t n=bs->GetVectorSize(), npts=R.size();
    mat_t<U> P(npts, n);
    const size_t nchunk=std::max<size_t>(1, npts/256);          // ~256 chunks: amortises, still balances
    const size_t nblk  =(npts+nchunk-1)/nchunk;
    auto fill=[&](size_t b)
    {
        const size_t lo=b*nchunk, hi=std::min(npts, lo+nchunk);
        if (lo>=hi) return;
        rvec3vec_t sub(hi-lo);
        for (size_t g=lo; g<hi; g++) sub[g-lo]=R[g];
        const mat_t<U> B=(*bs)(sub);
        for (size_t g=lo; g<hi; g++) for (size_t i=0; i<n; i++) P(g,i)=B(g-lo,i);
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1)
    {
        std::exception_ptr firstEx;                 // throw containment (an escape = std::terminate)
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t b=0; b<nblk; b++)
        {
            try { fill(b); }
            catch (...)
            {
                #pragma omp critical (xc_phi_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return P;
    }
#endif
    for (size_t b=0; b<nblk; b++) fill(b);
    return P;
}

// The two Collocation::Values overloads: same body, own cache per scalar.
template <class U> const mat_t<U>& DeltaFit_IBS::Table(std::map<Irrep,mat_t<U>>& cache, const Orbital_1E_IBS<U>& orb) const
{
    const Irrep id=orb.GetIrrep(Spin::None);       // SPATIAL key: the two spin channels share one table
    auto it=cache.find(id);
    if (it!=cache.end()) return it->second;
    mat_t<U> P=MakePhiAt<U>(&orb, itsQuad.mesh->Points(), "setup: XC-mesh Phi tables");
    ReportPhiSparsity(P, *itsQuad.mesh);           // mesh-aware (it reports SITE batching), hence outside the builder
    return cache.emplace(id, std::move(P)).first->second;
}

const mat_t<double>& DeltaFit_IBS::Values(const Orbital_1E_IBS<double>& orb) const {return Table(itsPhiR, orb);}
const mat_t<dcmplx>& DeltaFit_IBS::Values(const Orbital_1E_IBS<dcmplx>& orb) const {return Table(itsPhi , orb);}

} //namespace
