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
#include <memory>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <cmath>       // std::abs (the Phi-sparsity magnitude screen)
#include <complex>     // std::abs on a complex table entry
#include <map>
#include <exception>
#include <stdexcept>   // std::invalid_argument (the transition pair's shape checks)
#include <iostream>
#include <type_traits> // std::is_same_v (the real/complex contraction bodies)
#include <vector>
module qchem.BasisSet.DeltaFit_IBS;
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.Blaze;
import qchem.Reporting;   // report::Timed -- the setup timing bucket
import qchem.Parallel;    // WorkerThreads (QCHEM_OPENMP_THREADS): the table build is per-point-block parallel

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
    static const bool on=qchem::Diagnostics::Enabled("phi_sparsity");
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
    // question -- the table is bit-identical at any thread count.  Opt-in (QCHEM_OPENMP_THREADS) like every
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
    // and cannot host OpenMP, and molecular bases reach here too (MakeXCQuadrature serves every
    // structure), so handing the whole set to a serial default would deserialise them.  Each chunk writes its
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
    if (const int nthreads=qchem::InnerThreads(); nthreads>1)
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

// THE SHARED-Φ BATCHING (lane S 2026-11-19).  On the images realization a contraction forms Φ^k a point batch at a time from the shared
// record and GEMMs it; these two sizes shape that.  Neither is physics and neither moves a number beyond roundoff:
//  * kRowBatch: rows per batch.  Big enough that each batch's GEMM runs at BLAS speed on plain matrices (the whole-matrix-vs-view note at
//    AdjointDenseT: a plain destination keeps the zgemm), small enough that a batch's rows stay in cache (1024 x 152 x 16 B = 2.5 MB).
//  * kAdjointGroups: the adjoint's FIXED reduction partition.  Each group sums its batches in order into its own n x n partial, and the
//    partials are added in group order, so H is bit-identical at ANY thread count (the threads share out the groups, not the sum).
constexpr size_t kRowBatch=1024;
constexpr size_t kAdjointGroups=16;

// The shared image record of one lattice basis on my points: chunked and threaded exactly as the dense build (MakePhiAt), each chunk asked
// of the source in point order and appended in chunk order, so the record is identical at any thread count.  \a stride > 1 builds it on
// every stride-th point only: the policy's sample (ChoosePhiRealization).
static PointImages BuildImages(const BlochImageSource& src, const rvec3vec_t& R, size_t stride)
{
    std::vector<size_t> idx;
    for (size_t g=0; g<R.size(); g+=stride) idx.push_back(g);
    const size_t npts=idx.size();
    const size_t nchunk=std::max<size_t>(1, npts/256);
    const size_t nblk  =(npts+nchunk-1)/nchunk;      // every chunk non-empty
    std::vector<PointImages> parts(nblk);
    auto fill=[&](size_t b)
    {
        const size_t lo=b*nchunk, hi=std::min(npts, lo+nchunk);
        rvec3vec_t sub(hi-lo);
        for (size_t g=lo; g<hi; g++) sub[g-lo]=R[idx[g]];
        parts[b]=src.MakePointImages(sub);
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::InnerThreads(); nthreads>1)
    {
        std::exception_ptr firstEx;
        #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
        for (size_t b=0; b<nblk; b++)
        {
            try { fill(b); }
            catch (...)
            {
                #pragma omp critical (xc_images_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
    }
    else
#endif
        for (size_t b=0; b<nblk; b++) fill(b);
    if (nblk==0) throw std::logic_error("DeltaFit_IBS: an image record over an empty point set");
    PointImages all=std::move(parts[0]);
    for (size_t b=1; b<nblk; b++) {all.Append(parts[b]); parts[b]=PointImages();}
    return all;
}

template <class T> std::shared_ptr<const PointImages> DeltaFit_IBS<T>::Images(const BlochImageSource& src) const
{
    return itsImages.GetOrFill(src.ImageSource(), [this,&src]
    {
        qchem::report::Timed timed("setup: XC-mesh Phi images (one record, shared across k)");
        auto im=std::make_shared<const PointImages>(BuildImages(src, itsQuad.GetMesh()->Points(), 1));
        report::Accumulate("cache", "deltaPhiImageBytes", im->Bytes());   // the record reports itself as it fills, as the tables do
        return im;
    });
}

// The two Collocation::Values overloads: same body, own cache per scalar.  ONE OF TWO REALIZATIONS per block: the images (on the images
// route, for a block whose factory gave it an image source), else the dense table, filled exactly as before.
template <class T> template <class U> auto DeltaFit_IBS<T>::Table(const InstanceSlotMap<PhiRows<U>>& cache, const Orbital_1E_IBS<U>& orb) const
    -> const PhiRows<U>&
{
    // Thread-safe (2026-11-02, the threaded multi-k hang): every fitter and the density's forward reach these maps per block from the
    // parallel Fock pass; the slot map's once-fill builds each block's table once.
    return cache.GetOrFill(orb.Instance(), [this,&orb]
    {
        PhiRows<U> R;
        if (itsRealization==PhiRealization::Images && itsImageSourceOf)
            if (const BlochImageSource* src=itsImageSourceOf(orb); src && src->CanMakePointImages())
            {
                R.images=Images(*src);
                R.phases=src->ImagePhases(*R.images);
                if (R.images->NumFunctions()!=orb.GetNumFunctions() || R.images->NumPoints()!=itsQuad.GetMesh()->size())
                    throw std::logic_error("DeltaFit_IBS: a block's image record is not (my points x its functions)");
                report::Accumulate("cache", "deltaPhiImageBlocks", 1.0);
                return R;
            }
        mat_t<U> P;
        if constexpr (std::is_same_v<U,double>)   // the factory's own evaluation of its blocks (the atom's radial face), when it gave one
            if (itsRealValues)
            {
                qchem::report::Timed timed("setup: XC-mesh Phi tables");
                P=itsRealValues(orb, itsQuad.GetMesh()->Points());
                if (P.rows()!=itsQuad.GetMesh()->size() || P.columns()!=orb.GetNumFunctions())
                    throw std::logic_error("DeltaFit_IBS: the factory's point values for a block are not (points x functions)");
            }
        if (P.rows()==0) P=MakePhiAt<U>(&orb, itsQuad.GetMesh()->Points(), "setup: XC-mesh Phi tables");
        ReportPhiSparsity(P, *itsQuad.GetMesh());           // mesh-aware (it reports SITE batching), hence outside the builder
        report::Accumulate("cache", "deltaPhiBytes",   double(P.rows())*double(P.columns())*sizeof(U));   // S6: the table reports itself as it fills
        report::Accumulate("cache", "deltaPhiEntries", 1.0);
        R.dense=std::move(P);
        return R;
    });
}

template <class T> template <class U> mat_t<U> DeltaFit_IBS<T>::Rows(const PhiRows<U>& P, size_t lo, size_t hi) const
{
    if (P.IsDense()) return mat_t<U>(blazem::submatrix(P.dense, lo, size_t(0), hi-lo, P.dense.columns()));
    mat_t<dcmplx> C=P.images->Bloch(P.phases, lo, hi);
    if constexpr (std::is_same_v<U,dcmplx>) return C;
    else
    {
        mat_t<double> r(C.rows(), C.columns());     // a TRIM block: the phases are exactly ±1, so the imaginary parts are exactly 0 (Step 0)
        for (size_t i=0;i<C.rows();i++) for (size_t j=0;j<C.columns();j++) r(i,j)=std::real(C(i,j));
        return r;
    }
}

// THE FACTORY'S SETUP STEP (see the declaration).  The dense bytes are a count; the image bytes an exact count on a subsample of ~4096
// points (the image build costs one χ pass, about one dense block's build, so a full count would charge a Γ run a setup it has no use for),
// scaled to the whole mesh -- an ESTIMATE, labelled as one in the footprint.
template <class T> PhiRealization DeltaFit_IBS<T>::ChoosePhiRealization(const std::vector<const Orbital_1E_IBS<dcmplx>*>& cblocks,
                                                                        const std::vector<const Orbital_1E_IBS<double>*>&  rblocks,
                                                                        const PhiChooser& choose)
{
    const rvec3vec_t& pts=itsQuad.GetMesh()->Points();
    PhiFootprint f;
    f.points=pts.size();
    f.complexBlocks=cblocks.size();
    f.realBlocks=rblocks.size();
    bool all=bool(itsImageSourceOf);
    std::vector<const BlochImageSource*> sources;
    std::vector<uint64_t>                serials;    // a source's identity is its token's serial, never its address
    auto visit=[&](const IrrepBasisSet_IDs& b, size_t n, size_t scalarBytes)
    {
        f.denseBytes+=double(f.points)*double(n)*double(scalarBytes);
        f.maxFunctions=std::max(f.maxFunctions, n);
        const BlochImageSource* src=itsImageSourceOf ? itsImageSourceOf(b) : nullptr;
        if (!src || !src->CanMakePointImages()) {all=false; return;}
        const uint64_t id=src->ImageSource().Serial();
        if (std::find(serials.begin(), serials.end(), id)==serials.end()) {serials.push_back(id); sources.push_back(src);}
    };
    for (const auto* b : cblocks) visit(*b, b->GetNumFunctions(), sizeof(dcmplx));
    for (const auto* b : rblocks) visit(*b, b->GetNumFunctions(), sizeof(double));
    f.imagesAvailable=all && !sources.empty();
    f.imageSources=sources.size();
    if (f.imagesAvailable)
    {
        qchem::report::Timed timed("setup: XC-mesh Phi image count (subsample)");
        f.sampleStride=std::max<size_t>(1, f.points/4096);
        for (const auto* src : sources)
        {
            const PointImages smp=BuildImages(*src, pts, f.sampleStride);
            f.imagesBytes+=smp.Bytes()*double(f.points)/double(smp.NumPoints());
        }
        f.imagesBuildBytes=f.imagesBytes;     // the chunks, merged into one record, are freed as they are appended: at most one copy more
        const double n=double(f.maxFunctions), threads=double(std::max(1, qchem::InnerThreads()));
        const double nbat=std::ceil(double(f.points)/double(kRowBatch));
        f.imagesWorkBytes=threads*2.0*double(kRowBatch)*n*sizeof(dcmplx) + std::min(double(kAdjointGroups), nbat)*n*n*sizeof(dcmplx);
    }
    const PhiRealization r = choose ? choose(f) : PhiRealization::Dense;
    if (r==PhiRealization::Images && !f.imagesAvailable)
        throw std::logic_error("DeltaFit_IBS::ChoosePhiRealization: the chooser picked the images realization for blocks that have none");
    itsRealization=r;
    if (r==PhiRealization::Images)
        for (const auto* src : sources) Images(*src);     // built now, at setup (the run's console reads in run order)
    return r;
}

// THE 3-CENTRE OVERLAP <chi_i|delta_g|chi_j> = w_g conj(chi_i(r_g)) chi_j(r_g), as the house CONTRACTIBLE
// object: rank-3 and never materialised, so what the caller gets is the two contractions of it (three, once
// the density's factored form is counted -- one integral, two representations of D).  Phi and w stay in
// here; the caller supplies only what it owns, its coefficients or its density matrix.
template <class T> template <class U> const PairProjection<U>& DeltaFit_IBS<T>::Tensor(const InstanceSlotMap<PairProjection<U>>& cache,
                                                             const InstanceSlotMap<PhiRows<U>>& phis,
                                                             const Orbital_1E_IBS<U>& orb) const
{
    return cache.GetOrFill(orb.Instance(), [this,&phis,&orb]
    {
        // Bind the table ONCE here: a filled slot is address-stable, and the table and this record share the block's key, so they are
        // evicted together; the closures may capture it by reference and no contraction re-does the lookup.
        const PhiRows<U>& P=Table(phis, orb);
        PairProjection<U> g;
        g.applyRawAdjoint  =[this,&P](const rvec_t& c)   {return AdjointT<U>(P, c);};
        g.applyRaw         =[this,&P](const hmat_t<U>& D){return ForwardT<U>(P, D);};
        g.applyRawFactored =[this,&P](const mat_t<U>& L) {return ForwardFactoredT<U>(P, L);};
        return g;
    });
}

template <class T> const PairProjection<T>& DeltaFit_IBS<T>::Overlap3C(const Orbital_DFT_IBS<T,T>& orb) const
{
    if constexpr (std::is_same_v<T,double>) return Tensor(itsO3R, itsPhiR, orb);   // a block's table is typed by the BLOCK's scalar
    else                                    return Tensor(itsO3 , itsPhi , orb);
}

// THE TRANSITION PAIR (doc/LinearResponsePlan.md §3d B2, finding 4): the ground-state forward and adjoint with the
// bra's table on one side and the ket's on the other.  δρ_a = Σ_j [Φ^bra δD]_aj conj(Φ^ket_aj), and its adjoint
// Φ^bra† diag(w v) Φ^ket -- so Σ_ij conj(δD_ij) h_ij = Σ_a w_a conj(δρ_a) v_a exactly.  Complex point values (a
// Bloch-q function), not Hermitian, never folded.  Serial and unscreened: a response applies this a few hundred
// times per material, and the tables are the ground state's (efficiency is the deferred ledger's, §3d).
template <class T> cvec_t DeltaFit_IBS<T>::TransitionForwardImp(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                       const mat_t<dcmplx>& dD) const
{
    qchem::report::Timed timed("response: XC-mesh transition forward");
    const PhiRows<dcmplx>& Pb=Table(itsPhi, bra);
    const PhiRows<dcmplx>& Pk=Table(itsPhi, ket);
    if (dD.rows()!=Pb.NumCols() || dD.columns()!=Pk.NumCols())
        throw std::invalid_argument("DeltaFit_IBS::TransitionForward: δD is not bra x ket");
    const size_t npts=Pk.NumRows();
    cvec_t ro(npts);
    // Whole tables when both are dense (one GEMM, as before); a point batch at a time otherwise (the images realization).
    const size_t B = (Pb.IsDense() && Pk.IsDense()) ? std::max<size_t>(npts,1) : kRowBatch;
    for (size_t lo=0; lo<npts; lo+=B)
    {
        const size_t hi=std::min(npts, lo+B);
        const mat_t<dcmplx> PD = Pb.IsDense() && Pk.IsDense() ? mat_t<dcmplx>(Pb.dense*dD) : mat_t<dcmplx>(Rows(Pb,lo,hi)*dD);
        const mat_t<dcmplx> Kb = Pk.IsDense() ? mat_t<dcmplx>() : Rows(Pk,lo,hi);
        for (size_t g=lo; g<hi; g++)
        {
            dcmplx acc=0.0;
            for (size_t j=0; j<PD.columns(); j++) acc+=PD(g-lo,j)*std::conj(Pk.IsDense() ? Pk.dense(g,j) : Kb(g-lo,j));
            ro[g]=acc;
        }
    }
    return ro;
}
template <class T> mat_t<dcmplx> DeltaFit_IBS<T>::TransitionAdjointImp(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                              const cvec_t& v) const
{
    qchem::report::Timed timed("response: XC-mesh transition adjoint");
    const PhiRows<dcmplx>& Pb=Table(itsPhi, bra);
    const PhiRows<dcmplx>& Pk=Table(itsPhi, ket);
    const rvec_t& w=itsQuad.GetMesh()->Weights();
    const size_t npts=Pk.NumRows();
    if (v.size()!=npts) throw std::invalid_argument("DeltaFit_IBS::TransitionAdjoint: the field is not on my points");
    if (Pb.IsDense() && Pk.IsDense())      // the former single GEMM, unchanged
    {
        const mat_t<dcmplx>& K=Pk.dense;
        mat_t<dcmplx> WP(K.rows(), K.columns());
        for (size_t g=0; g<K.rows(); g++)
        {
            const dcmplx wv=w[g]*v[g];
            for (size_t j=0; j<K.columns(); j++) WP(g,j)=wv*K(g,j);
        }
        return mat_t<dcmplx>(blazem::trans(blazem::conj(Pb.dense))*WP);
    }
    mat_t<dcmplx> M(Pb.NumCols(), Pk.NumCols(), dcmplx(0.0));
    for (size_t lo=0; lo<npts; lo+=kRowBatch)       // a point batch at a time (the images realization), summed in batch order
    {
        const size_t hi=std::min(npts, lo+kRowBatch);
        const mat_t<dcmplx> Kb=Rows(Pk,lo,hi), Bb=Rows(Pb,lo,hi);
        mat_t<dcmplx> WP(hi-lo, Kb.columns());
        for (size_t g=lo; g<hi; g++)
        {
            const dcmplx wv=w[g]*v[g];
            for (size_t j=0; j<Kb.columns(); j++) WP(g-lo,j)=wv*Kb(g-lo,j);
        }
        M+=mat_t<dcmplx>(blazem::trans(blazem::conj(Bb))*WP);
    }
    return M;
}

// THE ADJOINT: <i|Sum_g c_g delta_g|j> = Phi^dag diag(w c) Phi -- scale the rows, one zgemm, hermitize.
// (The GEMM result is Hermitian up to roundoff; the explicit i<=j fill keeps chmat_t's invariant exactly.)
template <class T> template <class U> hmat_t<U> DeltaFit_IBS<T>::AdjointT(const PhiRows<U>& P, const rvec_t& v) const
{
    return P.IsDense() ? AdjointDenseT<U>(P.dense, v) : AdjointBatchedT<U>(P, v);
}

// The ADJOINT on the images realization: H = Σ_batches Φ_b† diag(w v)_b Φ_b, each batch's rows formed from the shared record, each batch's
// product a plain-matrix GEMM.  Threaded over the FIXED group partition (kAdjointGroups): bit-identical at any thread count.
template <class T> template <class U> hmat_t<U> DeltaFit_IBS<T>::AdjointBatchedT(const PhiRows<U>& P, const rvec_t& v) const
{
    const rvec_t& w=itsQuad.GetMesh()->Weights();
    const size_t npts=P.NumRows(), n=P.NumCols();
    if (v.size()!=npts) throw std::invalid_argument("DeltaFit_IBS::Adjoint: the field is not on my points");
    qchem::report::Timed timed("scf: XC-mesh quadrature H_xc (all iterations)");
    const size_t nbat=(npts+kRowBatch-1)/kRowBatch;
    const size_t ngrp=std::max<size_t>(1, std::min(kAdjointGroups, nbat));
    std::vector<mat_t<U>> part(ngrp, mat_t<U>(n, n, U(0.0)));
    auto group=[&](size_t gi)
    {
        const size_t b0=gi*nbat/ngrp, b1=(gi+1)*nbat/ngrp;
        for (size_t b=b0; b<b1; b++)
        {
            const size_t lo=b*kRowBatch, hi=std::min(npts, lo+kRowBatch);
            const mat_t<U> Pb=Rows(P, lo, hi);
            mat_t<U> WP(hi-lo, n);
            for (size_t q=0; q<hi-lo; q++)
            {
                const double wv=w[lo+q]*v[lo+q];
                for (size_t i=0; i<n; i++) WP(q,i)=wv*Pb(q,i);
            }
            part[gi]+=mat_t<U>(blazem::trans(blazem::conj(Pb))*WP);
        }
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::InnerThreads(); nthreads>1 && ngrp>1)
    {
        std::exception_ptr firstEx;
        #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
        for (size_t gi=0; gi<ngrp; gi++)
        {
            try { group(gi); }
            catch (...)
            {
                #pragma omp critical (xc_adj_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
    }
    else
#endif
        for (size_t gi=0; gi<ngrp; gi++) group(gi);
    mat_t<U> M=std::move(part[0]);
    for (size_t gi=1; gi<ngrp; gi++) M+=part[gi];
    hmat_t<U> H(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            H(i,j)=M(i,j);
    return H;
}

template <class T> template <class U> hmat_t<U> DeltaFit_IBS<T>::AdjointDenseT(const mat_t<U>& P, const rvec_t& v) const
{
    const rvec_t&        w=itsQuad.GetMesh()->Weights();  // my weights -- they never leave
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
    //    (qchem::FixBlasThreads, 1 by default): see that declaration for when raising it is safe -- at THIS
//    site there is no concurrent level above, so QCHEM_BLAS_THREADS>1 gives the zgemm the idle cores.
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
    if (const int nthreads=qchem::InnerThreads(); nthreads>1)
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
    const int nthreads=qchem::InnerThreads();
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
// functions.  Moved here VERBATIM from IrrepCD_Core::ProjectOnto (2026-08-23): it is a contraction of
// MY integral, so it belongs beside its adjoint, and having the two GEMMs in one file is also where the
// shared BLAS-dispatch reasoning above belongs.  Threaded by OUTPUT BLOCK -- row g depends on row g alone,
// so the result is bit-identical at any thread count.  Opt-in via QCHEM_OPENMP_THREADS (qchem.Parallel).
// The two forwards on either realization: the dense table keeps its own body (ForwardDenseT / the factored one below, unchanged); the
// images realization forms each point batch's rows and runs the same per-row reduction on them.  Row g depends on row g alone, so the
// batches thread with no reduction and the result is bit-identical at any thread count.
template <class T, class U, class F> static void ForEachRowBatch(size_t npts, F&& batch)
{
    const size_t nbat=(npts+kRowBatch-1)/kRowBatch;
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::InnerThreads(); nthreads>1 && nbat>1)
    {
        std::exception_ptr firstEx;
        #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
        for (size_t b=0; b<nbat; b++)
        {
            try { batch(b*kRowBatch, std::min(npts, (b+1)*kRowBatch)); }
            catch (...)
            {
                #pragma omp critical (xc_fwd_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return;
    }
#endif
    for (size_t b=0; b<nbat; b++) batch(b*kRowBatch, std::min(npts, (b+1)*kRowBatch));
}

template <class T> template <class U> rvec_t DeltaFit_IBS<T>::ForwardT(const PhiRows<U>& P, const hmat_t<U>& D) const
{
    if (P.IsDense()) return ForwardDenseT<U>(P.dense, D);
    if (P.NumCols()!=D.rows()) throw std::invalid_argument("DeltaFit_IBS::Forward: D is not over the block's functions");
    const size_t npts=P.NumRows();
    rvec_t ro(npts, 0.0);
    const mat_t<U> Dp(D);
    ForEachRowBatch<T,U>(npts, [&](size_t lo, size_t hi)
    {
        const mat_t<U> Pb=Rows(P, lo, hi);
        const mat_t<U> PD=Pb*Dp;
        for (size_t q=0; q<hi-lo; q++)
        {
            U acc{};
            for (size_t j=0; j<PD.columns(); j++)
                if constexpr (std::is_same_v<U,dcmplx>) acc+=PD(q,j)*std::conj(Pb(q,j));
                else                                    acc+=PD(q,j)*Pb(q,j);
            ro[lo+q]=std::real(acc);
        }
    });
    return ro;
}

template <class T> template <class U> rvec_t DeltaFit_IBS<T>::ForwardDenseT(const mat_t<U>& P, const hmat_t<U>& D) const
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
    if (const int nthreads=qchem::InnerThreads(); nthreads>1 && npts>size_t(nthreads))
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
template <class T> template <class U> rvec_t DeltaFit_IBS<T>::ForwardFactoredT(const PhiRows<U>& PR, const mat_t<U>& L) const
{
    if (!PR.IsDense())
    {
        if (PR.NumCols()!=L.rows()) throw std::invalid_argument("DeltaFit_IBS::ForwardFactored: L is not over the block's functions");
        const size_t npts=PR.NumRows();
        rvec_t ro(npts, 0.0);
        ForEachRowBatch<T,U>(npts, [&](size_t lo, size_t hi)
        {
            const mat_t<U> Psi=Rows(PR, lo, hi)*L;
            for (size_t q=0; q<hi-lo; q++)
            {
                double acc=0.0;
                for (size_t m=0; m<Psi.columns(); m++) acc+=std::norm(std::complex<double>(Psi(q,m)));
                ro[lo+q]=acc;
            }
        });
        return ro;
    }
    const mat_t<U>& P=PR.dense;
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
    if (const int nthreads=qchem::InnerThreads(); nthreads>1 && npts>size_t(nthreads))
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

// The two fit axes (2026-10-31): the periodic (Bloch) one and the molecular one.
template class DeltaFit_IBS<dcmplx>;
template class DeltaFit_IBS<double>;
// The real-TRIM record of the Bloch fit (RealTensor, inline in the interface, reaches it from other TUs): named here so it is emitted.
template const PairProjection<double>& DeltaFit_IBS<dcmplx>::Tensor<double>(const InstanceSlotMap<PairProjection<double>>&,
                                                                           const InstanceSlotMap<PhiRows<double>>&, const Orbital_1E_IBS<double>&) const;

} //namespace
