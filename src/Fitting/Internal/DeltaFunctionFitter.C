// File: Fitting/Internal/DeltaFunctionFitter.C  Scalar fitter on a DELTA (identity) fit basis.
//
// The third scalar fitter, and the one that makes the delta representation Liskov-substitutable with the
// other two: the XC term now obtains a fitter from the same Factory, calls the same DoFit/Overlap, and
// never learns which representation it got.
//
// WHAT IT IS.  A delta basis is ORTHOGONAL but not orthoNORMAL -- <delta_g|delta_h> = w_g delta_gh -- so
// the fit is the projection divided by the metric diagonal:
//
//     c_g = <delta_g|f> / <delta_g|delta_g> = w_g f(r_g) / w_g = f(r_g),
//
// i.e. the fit coefficients ARE the field's values at the mesh points -- but that is a FACT about this
// representation, not the shape of the code: the code is "projection, then my metric", the same two steps
// the Gaussian fitter takes with S^-1 and the orthoNORMAL sibling takes with nothing (the diagonal it
// asserts to be all ones -- OrthoNormalFunctionFitter.C).  Both halves are the BASIS's: it owns its
// functions and its metric, so nothing here knows where its points are or what kind of mesh they came from.
//
// AND WHY THE CONTRACTION IS THE BASIS'S.  Overlap is sum_g c_g <chi_i|delta_g|chi_j>, which is the
// weighted quadrature Phi^dag diag(w c) Phi -- the basis owns the points, the weights AND the Phi table,
// so it performs it and this class only decides WHICH c to hand over.  Both block scalars are served
// (FitContraction<double,dcmplx> and <dcmplx,dcmplx>), which a mixed real/complex run needs and which is
// exactly what the single-scalar fitter face could not express before the FitContraction split.
module;
#include <complex>   // std::real -- the Charge() narrow (2026-09-09)
#include <cassert>
#include <map>           // the per-block Phi handles (shallow Projector3 copies)
#include <memory>
#include <mutex>         // std::lock_guard over itsCacheLock (the per-block caches are filled from a parallel Fock pass)
#include <ostream>
#include <stdexcept>
#include <type_traits>   // std::is_same_v (same-scalar face vs the mixed-run cross-cast)
export module qchem.Fitting.Internal.DeltaFunctionFitter;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar + FitContraction<U,TFit>, ProjectedScalar_R
import qchem.Fitting.Types;                   // robs_t<U>
import qchem.BasisSet.Orbital_DFT_IBS;                // cFIT_SF_ABS + Integrals_Overlap3C<U,TFit> (the held fit basis)
import qchem.Fitting.FitOperations;           // OrthogonalFit -- this library's projection/diagonal fit
import qchem.Blaze;                           // hmat_t<U>
import qchem.Parallel;                        // LazyLock -- the per-block caches' guard

export namespace qchem::Fitting
{

//! \brief THE PERIODIC FIT AXIS'S EXTRA FACES (2026-10-31), which only a fitter over a COMPLEX (Bloch) \f$\delta\f$ basis has:
//! the real TRIM block's contraction (3c-3, \c FitContraction<double,dcmplx>) and the FORWARD vendor the periodic density
//! projects onto (\c ScalarProjector, whose faces are the Bloch blocks').  A base that exists only for \c dcmplx (CRTP: the
//! bodies are the fitter's own), so the real (molecular) instantiation carries no face it cannot serve.
template <class Self, class T> class DeltaPeriodicFitFaces {};
template <class Self> class DeltaPeriodicFitFaces<Self,dcmplx>
    : public virtual ScalarProjector                 // ...the FORWARD direction, off the same handles
    , public virtual FitContraction<double,dcmplx>   // a real TRIM block's quadrature runs in real arithmetic
{
public:
    virtual hmat_t<double> Overlap(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return Self_().template Contract<double>(orb);}
    //! \copydoc Fitting::ScalarProjector::Forward
    //! Built over the SAME held handle the adjoint contracts -- which is the point of holding it: the density's forward and
    //! the term's adjoint come off ONE constructed object.  ▶ The fitter is "ACTOR 2": it builds that object per block and
    //! hands each client only its own half -- `IrrepCD` the FORWARD (through this face), the XC term the ADJOINT (\c Overlap).
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override
        {return Self_().template Integrator<double>(orb);}
    const qcMesh::MatrixForward<dcmplx>& Forward(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override
        {return Self_().template Integrator<dcmplx>(orb);}
private:
    const Self& Self_() const {return static_cast<const Self&>(*this);}
};

//! \brief Scalar (overlap-metric) fitter on a \f$\delta\f$ fit basis -- the projection IS the fit, up to
//! the basis's own metric diagonal.  ONE class template on the fit axis's scalar \a T (2026-10-31, as \c DeltaFit_IBS<T>):
//! \c dcmplx on the periodic path (serving BOTH orbital-block scalars, as a mixed run requires -- the extra faces come from
//! \c DeltaPeriodicFitFaces), \c double on the molecular one.
template <class T> class DeltaScalarFitter
    : public virtual FunctionFitter_Scalar
    , public virtual FitContraction<T,T>             // the same-scalar blocks' contraction
    , public virtual FieldProjector                  // the field projection, on either fit axis
    , public         DeltaPeriodicFitFaces<DeltaScalarFitter<T>,T>   // dcmplx only: the real TRIM block + the forward vendor
{
    friend class DeltaPeriodicFitFaces<DeltaScalarFitter<T>,T>;
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS<T>> fbs_t;
    explicit DeltaScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs)
    {
        assert(itsFitBasis && itsFitBasis->isOrtho() &&
               "DeltaScalarFitter: a delta basis is orthogonal (diagonal metric) by construction");
    }

    //! The fit, in the ONE shape every scalar fitter has: take the basis's projection and apply my metric.
    //! Here the metric is diagonal, so \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$ -- which for
    //! \f$\delta\f$ is \f$w_g f(r_g)/w_g\f$.  \c OrthogonalFit is that one line, shared with the
    //! only other consumer of it (a matrix-free density expressing itself over a \f$\delta\f$ basis).
    //! \note This is where \c OverlapDiagonal() finally gets a production consumer.  The division does NOT
    //! cancel to the last bit -- \f$\mathrm{fl}(w_gf_g)/w_g\ne f_g\f$ in general -- so the coefficients can
    //! move by an ulp against the old direct sample.  Measured on the pinned Si two-route gate: no printed
    //! digit and no iteration count moved (doc/CleanupCandidates.md R1.0 increment 3).
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        itsC=OrthogonalFit(*itsFitBasis, *ps.GetScalarFunction());
    }
    //! \copydoc Fitting::FieldProjector::Project
    //! Structurally identical to \c DoFit and deliberately NOT sharing its buffer: the same projection, returned.
    rvec_t Project(const ScalarFunction<double>& f) const override {return OrthogonalFit(*itsFitBasis, f);}
    //! \copydoc Fitting::FieldProjector::NumCoefficients
    size_t NumCoefficients() const override {return itsFitBasis->GetNumFunctions();}

    //! \copydoc Fitting::FitContraction::Overlap
    //! \f$\langle\chi_i|\sum_g c_g\delta_g|\chi_j\rangle\f$: ask the basis for the 3-centre overlap over
    //! ITS functions and contract it with MY coefficients.  Structurally the same two calls the Gaussian
    //! fitter makes (\c Overlap3C then contract) -- which is what Liskov substitutability of the
    //! representations means in practice, rather than by analogy.
    virtual hmat_t<T> Overlap(const BasisSet::Orbital_DFT_IBS<T,T>& orb) const override {return Contract<T>(orb);}

    virtual void ReScale(double factor) override {itsC *= factor;}

    //! \warning NO evaluatable fitted field, by design.  A \f$\delta\f$ expansion has no value BETWEEN its
    //! points -- \f$\sum_g c_g\delta_g(r)\f$ is a distribution -- so this fitter does not derive the
    //! \c ScalarFunction capability its two siblings do, and nothing asks it to.  (The AO fitter evaluates
    //! its Gaussians, the plane-wave fitter inverse-transforms its coefficients; both are genuine fields.)
    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "DeltaScalarFitter (delta basis: the fit coefficients ARE the sampled values)" << std::endl;}

private:
    //! ONE body for both block scalars: the basis's \f$\langle\chi_i|\delta_g|\chi_j\rangle\f$, contracted
    //! against my coefficients.  A mixed run (3c-3) reaches it with either scalar.
    //! ★ THROUGH THE INTEGRATOR, NOT THE TENSOR (2026-09-09) -- the ADJOINT half of the R1.0o migration,
    //! which had stopped after the forward.
    //!
    //! This used to be `Handle(...)` then `o3.applyRawAdjoint(itsC)` behind an `assert(o3.applyRawAdjoint)`:
    //! reaching into a raw `std::function` on the CONCRETE tensor, guarded by a check that is compiled out
    //! under NDEBUG -- i.e. in every production run.  Exactly the pattern `Forward` above stopped doing,
    //! five lines away, on the SAME held handle.  ⇒ Both directions now come off ONE
    //! object, so the sentence in `Forward`'s comment ("the density's forward and the term's adjoint
    //! provably come off ONE tensor") is finally something the CODE enforces rather than something the
    //! comment asserts.  The absent-pair check moved with it, from a dead assert to the
    //! integrator's constructor, which THROWS.
    //!
    //! The 3-centre overlap is a DFT-TIER question, and since 2026-08-24 the CONTRACTION face says so
    //! too: both overloads above take the DFT block, so the cross-cast that used to stand here belongs
    //! to (and now lives with) the caller.  TFit is \a T, the fit axis's scalar; on the periodic axis a real TRIM block
    //! against it is exactly the 3c-3 mixed case (U = double, T = dcmplx).
    template <class U> hmat_t<U> Contract(const BasisSet::Orbital_DFT_IBS<U,T>& orb) const
    {
        const qcMesh::MatrixAdjoint<U>& adj=Integrator<U>(orb);
        return adj.Adjoint(itsC);
    }
    //! \brief THE SHALLOW Φ HANDLE for one orbital block, held for the run (user, 2026-08-24).
    //!
    //! \c Projector3 is three closures over the fit basis's OWN address-stable \f$\Phi\f$ table, so a copy
    //! is ~100 bytes and the npts×n table is never duplicated -- which is what makes "hold it" affordable
    //! at all.  Held rather than re-fetched so BOTH directions demonstrably come off one object.
    //! \warning the closures capture the producing basis: valid only while \c itsFitBasis lives, which it
    //! does, because this class co-owns it.
    //! \pre the caller holds \c itsCacheLock (\c Integrator, its only caller).
    template <class U> const Projector3<U>& Handle(SymMap<Projector3<U>>& cache,
                                                   const BasisSet::Orbital_DFT_IBS<U,T>& orb) const
    {
        const sym_t& id=orb.GetSymt();             // the block's SPATIAL symmetry, as the basis's own table cache uses
        auto it=cache.find(id);
        if (it!=cache.end()) return it->second;
        // The genuine 3c-3 "can you serve a REAL block?" question for U==double; for U==dcmplx it is a
        // derived-to-virtual-base cast that cannot fail.  Asked ONCE per block per run, here.
        const auto& face=dynamic_cast<const BasisSet::Integrals_Overlap3C<U,T>&>(*itsFitBasis);
        return cache.emplace(id, face.Overlap3C(orb)).first->second;   // shallow: closures, not the table
    }
    //! Which typed cache \a U uses -- so \c Contract and \c Overlap3C name one place, not two.
    template <class U> SymMap<Projector3<U>>& Cache() const
    {
        if constexpr (std::is_same_v<U,double>) return itsO3R;
        else                                    return itsO3;
    }
    //! Which typed integrator cache \a U uses -- the sibling of \c Cache, so the two never drift apart.
    template <class U> SymMap<ScreenedMatrixIntegrator<U>>& IntCache() const
    {
        if constexpr (std::is_same_v<U,double>) return itsIntR;
        else                                    return itsInt;
    }
    //! \brief The forward+adjoint VIEW of each held handle, cached beside it and keyed identically -- BOTH
    //! directions, since 2026-09-09 (it vended only the forward before; see \c Contract).  Stored BY VALUE
    //! in the map rather than as an owning pointer: the lifetime question then has no answer to get wrong,
    //! which is what CLAUDE.md's ownership rule is for.
    //!
    //! It must be CACHED, not returned by value: the faces hand back a reference, and a
    //! \c ScreenedMatrixIntegrator borrows its \c Projector3, so a temporary would dangle the moment the
    //! caller used it.  Keyed on the same spatial \c Symmetry as the tensor, and created only after the tensor is,
    //! so the two maps never disagree about which blocks exist.
    //!
    //! ⛔ THREAD-SAFE BY A LOCK, NOT BY A WARM-UP (2026-11-02, OpenWork §3 "THREADED multi-k run HANGS FOREVER").  The XC sampler
    //! builds ONE fitter per field value v (\c tSinglesDensitySampler::FitFor) and every Bloch block of the parallel Fock pass
    //! contracts against it, so these maps are FIRST filled inside that pass -- no serial cold assembly can warm a fitter that
    //! did not exist yet.  Unlocked, two blocks' find/emplace raced and corrupted the tree: bcc Li k8 at 8 threads spun in
    //! std::map::find for 10 h.  The lock is per call (one per block per iteration), and a reference into a std::map survives
    //! later insertions, so what is returned stays valid after the lock is released.
    template <class U> const ScreenedMatrixIntegrator<U>&
    Integrator(const BasisSet::Orbital_DFT_IBS<U,T>& orb) const
    {
        std::lock_guard<std::mutex> lock(itsCacheLock.m);   // guards itsO3R/itsO3/itsIntR/itsInt: Handle is reached only from here
        auto& fwd=IntCache<U>();
        const sym_t& id=orb.GetSymt();
        auto it=fwd.find(id);
        if (it!=fwd.end()) return it->second;
        const Projector3<U>& g=Handle(Cache<U>(), orb);   // materialise (or find) the tensor first
        // The integrator's "weights" are the FUNCTION integrals <f_a|1> of my coefficient axis, not point
        // weights -- which on a delta basis are numerically the mesh weights, and on any other fit basis
        // are still the right thing.  Charge() is that vector; it is real for a real-valued fit basis, and
        // the narrow is where that fact is stated once.
        const vec_t<T> q=itsFitBasis->Charge();
        rvec_t w(q.size());
        for (size_t a=0;a<q.size();a++) w[a]=std::real(q[a]);
        return fwd.emplace(id, ScreenedMatrixIntegrator<U>(g, std::move(w))).first->second;
    }

    mutable SymMap<Projector3<double>> itsO3R;   //!< real TRIM blocks' handles (3c-3), by SPATIAL symmetry
    mutable SymMap<Projector3<dcmplx>> itsO3;    //!< Bloch blocks' handles
    //! Their forward+adjoint views -- BOTH halves off one object per block (2026-09-09).
    mutable SymMap<ScreenedMatrixIntegrator<double>> itsIntR;
    mutable SymMap<ScreenedMatrixIntegrator<dcmplx>> itsInt;
    qchem::LazyLock itsCacheLock;   //!< the four maps above: filled from the parallel Fock pass (see Integrator)
    fbs_t  itsFitBasis;   //!< the δ basis -- my functions, their metric, and their 3-centre overlap
    rvec_t itsC;          //!< MY fit coefficients over that basis (see DoFit)
};

} //namespace
