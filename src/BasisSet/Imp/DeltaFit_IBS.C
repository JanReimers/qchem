// File: BasisSet/Imp/DeltaFit_IBS.C  The delta basis's collocation tables.
//
// THE 3-CENTRE OVERLAP <chi_i|delta_g|chi_j> = w_g conj(chi_i(r_g)) chi_j(r_g) -- job (2) of a fit basis
// (user, 2026-08-23): provide the integrals over ITS OWN functions that a Hamiltonian term needs to form
// H_ij and E.  Both contraction directions live here, beside each other:
//
//   ADJOINT  <i|Sum_g c_g delta_g|j>   -- the term/fitter's side; c comes from the FITTER, which is the
//                                         object that holds a fit
//   FORWARD  <delta_g|rho[D]>/w_g      -- the density's side; D comes from the DENSITY, in whichever of
//                                         its two forms (full, or a thin factor)
//
// Phi_gi = chi_i(r_g) is HOW this class evaluates that integral -- in principle op(r) from the orbital
// basis, in practice a cached table (user) -- and is entirely private: a table of values is not an
// integral and has no business in an interface.  Cached per BLOCK and geometry-fixed: built once per run
// per block, serving every SCF iteration.
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
#include <type_traits> // std::is_same_v (the real/complex contraction bodies)
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

// THE 3-CENTRE OVERLAP <chi_i|delta_g|chi_j> = w_g conj(chi_i(r_g)) chi_j(r_g), as the house CONTRACTIBLE
// object: rank-3 and never materialised, so what the caller gets is the two contractions of it (three, once
// the density's factored form is counted -- one integral, two representations of D).  Phi and w stay in
// here; the caller supplies only what it owns, its coefficients or its density matrix.
template <class U> const Projector3<U>& DeltaFit_IBS::Tensor(std::map<Irrep,Projector3<U>>& cache,
                                                             std::map<Irrep,mat_t<U>>& phis,
                                                             const Orbital_1E_IBS<U>& orb) const
{
    const Irrep id=orb.GetIrrep(Spin::None);       // SPATIAL key, as for the table itself
    auto it=cache.find(id);
    if (it!=cache.end()) return it->second;
    // Bind the table ONCE here: std::map nodes are address-stable, so the closures may capture it by
    // reference and no contraction re-does the lookup.
    const mat_t<U>& P=Table(phis, orb);
    Projector3<U> g;
    g.applyRawAdjoint  =[this,&P](const rvec_t& c)   {return AdjointT<U>(P, c);};
    g.applyRaw         =[this,&P](const hmat_t<U>& D){return ForwardT<U>(P, D);};
    g.applyRawFactored =[this,&P](const mat_t<U>& L) {return ForwardFactoredT<U>(P, L);};
    return cache.emplace(id, std::move(g)).first->second;
}

const Projector3<double>& DeltaFit_IBS::Overlap3C(const Orbital_DFT_IBS<double,dcmplx>& orb) const {return Tensor(itsO3R, itsPhiR, orb);}
const Projector3<dcmplx>& DeltaFit_IBS::Overlap3C(const Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const {return Tensor(itsO3 , itsPhi , orb);}

// THE ADJOINT: <i|Sum_g c_g delta_g|j> = Phi^dag diag(w c) Phi -- scale the rows, one zgemm, hermitize.
// (The GEMM result is Hermitian up to roundoff; the explicit i<=j fill keeps chmat_t's invariant exactly.)
template <class U> hmat_t<U> DeltaFit_IBS::AdjointT(const mat_t<U>& P, const rvec_t& v) const
{
    const rvec_t&        w=itsQuad.mesh->Weights();  // my weights -- they never leave
    assert(v.size()==P.rows());
    qchem::report::Timed timed("scf: XC-mesh quadrature H_xc (all iterations)");
    mat_t<U> WP(P.rows(), P.columns());
    auto scale=[&](size_t g)
    {
        const double wv=w[g]*v[g];
        for (size_t i=0; i<P.columns(); i++) WP(g,i)=wv*P(g,i);
    };
    const size_t n=P.columns(), npts=P.rows();
    mat_t<U> M(n, n, U(0.0));
    // The per-iteration quadrature GEMM (once per XC term per spin): O(npts n^2), and with the rho
    // sampling dispatched to BLAS this is the SCF loop's largest bucket.  The row scaling above is
    // trivially parallel (row g is private to g) either way; the PRODUCT has two shapes:
    //
    //  * WITH BLAZE_BLAS_MODE -- ONE whole-matrix product, NO blocking and NO OpenMP.  Blaze hands a
    //    product to zgemm only when the DESTINATION is a plain matrix: assigning into a submatrix view
    //    silently drops it back to Blaze's own kernel.  Measured on this exact shape (48033x122,
    //    1 thread): whole-matrix 34.1 GFlop/s vs 1.87 for ANY blocked/viewed form -- so blocking to
    //    halve the flops (Hermitian) or to spread over 8 threads LOSES 13x to save 2x or 8x.  One
    //    dispatched zgemm beats every hand-parallel arrangement here, and it composes with the pin
    //    (qchem::PinBlasToOneThread): our parallelism stays at the levels above.
    //  * WITHOUT it -- the TRIANGULAR, output-blocked, threaded form: H is Hermitian and chmat_t
    //    stores i<=j, so column block [j0,j1) needs only rows [0,j1) of the left operand.  Every
    //    element M(i,j) is still ONE dot product over ALL mesh points accumulated by ONE thread in the
    //    serial order (a partition of the OUTPUT, not of the reduction), so it is bit-identical at any
    //    thread count.  Blocks carry unequal work, hence dynamic scheduling over ~4 blocks per thread.
    //
    // Either way only the upper triangle is READ below; the retired `½(M+M†)` symmetrisation acted on
    // roundoff alone (w and v are real, so M is Hermitian in exact arithmetic).
#ifdef BLAZE_BLAS_MODE
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1)
    {
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t g=0; g<npts; g++) scale(g);
    }
    else
#endif
        for (size_t g=0; g<npts; g++) scale(g);
    M = blazem::trans(blazem::conj(P))*WP;                    // one zgemm; see the note above
#else
    auto triBlock=[&](size_t j0, size_t nj)
    {
        const size_t j1=j0+nj;
        blazem::submatrix(M,0,j0,j1,nj) = blazem::trans(blazem::conj(blazem::submatrix(P,0,size_t(0),npts,j1)))
                                        * blazem::submatrix(WP,0,j0,npts,nj);
    };
#ifdef QCHEM_OPENMP
    const int nthreads=qchem::WorkerThreads();
    if (nthreads>1)
    {
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t g=0; g<npts; g++) scale(g);
        const size_t blk=std::max<size_t>(1, (n+4*size_t(nthreads)-1)/(4*size_t(nthreads)));
        const size_t nb=(n+blk-1)/blk;
        #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
        for (size_t b=0; b<nb; b++) triBlock(b*blk, std::min(blk, n-b*blk));
    }
    else
#endif
    {
        for (size_t g=0; g<npts; g++) scale(g);
        triBlock(0, n);
    }
#endif
    hmat_t<U> H(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            H(i,j)=M(i,j);
    return H;
}

// THE FORWARD: <delta_g|rho[D]>/w_g = [Phi D Phi^dag]_gg -- rho's expansion coefficients over my
// functions.  Moved here VERBATIM from IrrepCD_Core::DM_RhoAtPoints (2026-08-23): it is a contraction of
// MY integral, so it belongs beside its adjoint, and having the two GEMMs in one file is also where the
// shared BLAS-dispatch reasoning above belongs.  Threaded by OUTPUT BLOCK -- row g depends on row g alone,
// so the result is bit-identical at any thread count.  Opt-in via GPW_OMP_THREADS (qchem.Parallel).
template <class U> rvec_t DeltaFit_IBS::ForwardT(const mat_t<U>& P, const hmat_t<U>& D) const
{
    const size_t npts=P.rows();
    assert(P.columns()==D.rows());
    rvec_t ro(npts, 0.0);
    // D as a PLAIN matrix: an adaptor operand (HermitianMatrix) is not BLAS-dispatchable, so under
    // BLAZE_BLAS_MODE the product would silently fall back to blaze's own kernel -- 16x slower than zgemm
    // on this shape (measured).  One n x n copy per call is nothing beside the npts x n x n GEMM.
    const mat_t<U> Dp(D);
    auto block=[&](size_t g0, size_t g1)
    {
        auto Pb = blazem::submatrix(P, g0, size_t(0), g1-g0, P.columns());
        mat_t<U> PD = Pb*Dp;
        for (size_t g=g0; g<g1; g++)
        {
            U acc{};
            for (size_t j=0; j<PD.columns(); j++)
                if constexpr (std::is_same_v<U,dcmplx>) acc+=PD(g-g0,j)*std::conj(P(g,j));
                else                                    acc+=PD(g-g0,j)*P(g,j);
            ro[g]=std::real(acc);
        }
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1 && npts>size_t(nthreads))
    {
        const size_t blk=(npts+size_t(nthreads)-1)/size_t(nthreads);
        std::exception_ptr firstEx;                   // throw containment (an escape = std::terminate)
        #pragma omp parallel for schedule(static,1) num_threads(nthreads)
        for (size_t b=0; b<(npts+blk-1)/blk; b++)
        {
            try { block(b*blk, std::min(npts, (b+1)*blk)); }
            catch (...)
            {
                #pragma omp critical (cd_rho_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return ro;
    }
#endif
    // Serial: multiply the WHOLE table, not a full-height view -- a submatrix operand costs blaze its
    // aligned/padded kernel choice, and this is the default (unthreaded) path every suite run takes.
    mat_t<U> PD = P*Dp;
    for (size_t g=0; g<npts; g++)
    {
        U acc{};
        for (size_t j=0; j<PD.columns(); j++)
            if constexpr (std::is_same_v<U,dcmplx>) acc+=PD(g,j)*std::conj(P(g,j));
            else                                    acc+=PD(g,j)*P(g,j);
        ro[g]=std::real(acc);
    }
    return ro;
}

// THE SAME FORWARD for a caller holding D in FACTORED form, D = L L^dag: rho_g = ||[Phi L]_g||^2, a thin
// (npts x r) GEMM plus a row norm instead of the square one.  Non-negative BY CONSTRUCTION.  WHICH form
// the caller has -- and whether the rank pays -- is its business (ChargeDensity::FactoredRho owns the
// factor, its memo and the pays-for-itself test); this only knows how to contract against my table.
template <class U> rvec_t DeltaFit_IBS::ForwardFactoredT(const mat_t<U>& P, const mat_t<U>& L) const
{
    const size_t npts=P.rows();
    assert(P.columns()==L.rows());
    rvec_t ro(npts, 0.0);
    auto block=[&](size_t g0, size_t g1)
    {
        auto Pb = blazem::submatrix(P, g0, size_t(0), g1-g0, P.columns());
        mat_t<U> Psi = Pb*L;                          // (npts_b x rank) -- the thin GEMM
        for (size_t g=g0; g<g1; g++)
        {
            double acc=0.0;
            for (size_t m=0; m<Psi.columns(); m++) acc+=std::norm(std::complex<double>(Psi(g-g0,m)));
            ro[g]=acc;
        }
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1 && npts>size_t(nthreads))
    {
        const size_t blk=(npts+size_t(nthreads)-1)/size_t(nthreads);
        std::exception_ptr firstEx;
        #pragma omp parallel for schedule(static,1) num_threads(nthreads)
        for (size_t b=0; b<(npts+blk-1)/blk; b++)
        {
            try { block(b*blk, std::min(npts, (b+1)*blk)); }
            catch (...)
            {
                #pragma omp critical (cd_rho_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return ro;
    }
#endif
    block(0, npts);
    return ro;
}

} //namespace
