// File: HamiltonianImp.C  General matrix implementation of a Hamiltonian operator.
module;
#include <limits>
#include <vector>
#include <memory>
#include <string>
#include <type_traits>
export module qchem.Hamiltonian.Internal.Hamiltonian;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Hamiltonian.Internal.RealSpaceResidual;
import qchem.Structure;

export namespace qchem::Hamiltonian
{

//! \brief The REAL-BLOCK assembly (Step 3c-2) of a COMPLEX-run Hamiltonian: fold the terms' real block faces
//! (tStatic_CC<double> / tDynamic_CC<double>) into one real Fock/KS block for a real TRIM basis block -- the
//! \c Ham_RealBlock face the real WF children drive.  CRTP on the implementation (\a Imp, whose term lists it reads),
//! the shape of \c Composite_HFSystem (ChargeDensity/CompositeCD.C): the face exists on \c tHamiltonian<dcmplx> only,
//! so its implementation is a base of \c tHamiltonianImp<dcmplx> only, and the real Hamiltonian declares nothing it
//! cannot override (it once carried this member unconditionally, overriding nothing on <double>).
template <class Imp> class HamiltonianImp_RealBlock : public virtual tHamiltonian<dcmplx>
{
public:
    virtual hmat_t<double> GetMatrix(const tobs_t<double>*,const Spin& S,const tChargeDensity<dcmplx>*,
                                     const tbs_t<dcmplx>* wholeBasis) const override;
private:
    const Imp& self() const {return static_cast<const Imp&>(*this);}
};
struct NoHamiltonianImp_RealBlock {};
template <class T, class Imp> using HamiltonianImpRealBlockBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, HamiltonianImp_RealBlock<Imp>, NoHamiltonianImp_RealBlock>;

// Templated on T (rX/cX convention): the bare rHamiltonianImp is the <double> alias the existing concrete
// Hamiltonians derive from; cHamiltonianImp is the dcmplx (plane-wave) instantiation.
template <class T> class tHamiltonianImp
    : public virtual tHamiltonian<T>
    , public HamiltonianImpRealBlockBase<T,tHamiltonianImp<T>>   // the real-block assembly: complex instantiation only
{
    friend class HamiltonianImp_RealBlock<tHamiltonianImp<T>>;   // it folds this Hamiltonian's own term lists
public:
    explicit tHamiltonianImp(SpinGroup g);   //!< the imposed spin subgroup this Hamiltonian is built for
    virtual void Add(   tStatic_HT<T>* ) override;
    virtual void Add(  tDynamic_HT<T>*) override;
    virtual void Add(tDynamic_HF_HT<T>*) override;

    using tHamiltonian<T>::GetMatrix;   // keep the base 3-arg (null-basis) convenience overload visible
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>*,const Spin& S,const tChargeDensity<T>*,const tbs_t<T>* wholeBasis) const override;
    //! The REAL-BLOCK assembly lives in \c HamiltonianImp_RealBlock (a base on the complex instantiation only).
    virtual EnergyBreakdown GetTotalEnergy  (const tDM_CD<T>* ) const override;
    virtual std::vector<TermEnergy>     MakeTermEnergies(const tDM_CD<T>*) const override;
    virtual std::vector<hmat_t<T>>      MakeTermMatrices(const tobs_t<T>*, const Spin&, const tChargeDensity<T>*) const override;
    virtual std::vector<hmat_t<double>> MakeTermMatrices(const tobs_t<double>*, const Spin&, const tChargeDensity<dcmplx>*) const override;
    virtual SpinGroup       GetSpinGroup() const override {return itsSpinGroup;}
    virtual bool            IsRelativistic() const override {return itsIsRelativistic;}
    //! CONJUNCTIVE over the terms -- see the base declaration.  One non-Coulombic term (a PP) invalidates
    //! the virial for the whole Hamiltonian, so this is AND-ed in Add() where the two flags above are OR-ed.
    virtual bool            IsVirialValid () const override {return itsIsVirialValid;}
    //! CONJUNCTIVE like IsVirialValid: one SOC / vector-potential term flips every block complex
    //! (doc/RealComplexPlan.md); AND-ed in Add().
    virtual bool            PreservesReal () const override {return itsPreservesReal;}
    //! A density MATRIX (not a matrix-free fit) is required to seed the SCF iff this Hamiltonian either holds
    //! an exact-exchange HF term (K needs D) or is relativistic (the non-rel LDA-sibling SAD bootstrap does
    //! not apply -- see SCFIterator).  Derived from the term lists, so no concrete Hamiltonian need declare it.
    virtual bool            RequiresDensityMatrix() const override {return !itsHF_HTs.empty() || itsIsRelativistic;}
    //! FIRST term that owns an atom-centred partition wins (there is exactly one XC quadrature in a term
    //! stack; a second would be a construction error, not a case to average).
    virtual rvec_t          SiteMoments(const tChargeDensity<T>* cd) const override
    {
        for (const auto& t : itsDHTs)
            if (rvec_t m=t->SiteMoments(cd); m.size()>0) return m;
        return rvec_t();
    }
    //! \copydoc tHamiltonian::DeltaRho
    virtual double          DeltaRho(const tChargeDensity<T>* in, const tChargeDensity<T>* out) const override;
    //! \copydoc tHamiltonian::DeltaRhoLabel
    virtual std::string     DeltaRhoLabel() const override {return itsRealSpaceTerm ? "Δρ" : "‖ΔD‖";}
    //! \copydoc tHamiltonian::MakeACBN0Estimator
    //! The FIRST dynamic term answering the \c HubbardProjection face (abstract->abstract); null when none does.
    virtual std::unique_ptr<ACBN0Estimator> MakeACBN0Estimator() const override;
    //! \copydoc tHamiltonian::GetHubbardChannels
    virtual const HubbardChannels* GetHubbardChannels() const override;
    //! \copydoc tHamiltonian::GetHubbardUTarget
    virtual HubbardUTarget* GetHubbardUTarget() override;
    virtual std::unique_ptr<ResponseKernel<T>> MakeResponseKernel(const tbs_t<T>* wholeBasis,
                                                                  const tChargeDensity<T>* D0) const override;
    // The force fold (GeometryOptimizationPlan §13): Imp/HamiltonianForces.C.
    virtual std::vector<std::string> MissingGradients(const tbs_t<T>* wholeBasis, const Structure* st) const override;
    virtual rvec3vec_t MakeForces(const tbs_t<T>* wholeBasis, const tDM_CD<T>* D, const tDM_CD<T>* W, const Structure* st) const override;
    virtual std::vector<std::string> MissingVirials(const tbs_t<T>* wholeBasis, const Structure* st) const override;
    virtual rmat3d_t MakeVirial(const tbs_t<T>* wholeBasis, const tDM_CD<T>* D, const tDM_CD<T>* W, const Structure* st) const override;
    virtual std::vector<TermForce> MakeTermForces(const tbs_t<T>* wholeBasis, const tDM_CD<T>* D, const Structure* st) const override;
    virtual std::vector<TermVirial> MakeTermVirials(const tbs_t<T>* wholeBasis, const tDM_CD<T>* D, const Structure* st) const override;
    //! THE EAGER REFRESH PHASE (doc/OpenWork.md item **KP**): pre-warm every term's k-independent,
    //! density-dependent memo, once, before the caller's per-block Fock loop.  Folded over the DYNAMIC
    //! terms only -- a static term is density-independent by definition, so it has nothing to warm.
    //! Priced in its own report bucket by the caller, not here (this face is I/O-free).
    //! ⚠ TWO DUTIES SINCE 2026-09-09 (R1.0h), and they are one phase: everything that must happen before
    //! the block loop so the loop can be read-only.
    //!  1. PRE-CREATE this iteration's per-irrep cache slots (needs \a bs -- the block list).
    //!  2. PRE-WARM the k-independent density-derived memos (needs \a cd).
    //! Duty 1 runs over STATIC terms too, and that does not breach the "the phase must never reach a static
    //! term" rule: that rule is about refreshing FOR A DENSITY, which a static term has no business doing.
    //! Creating its slots is density-INDEPENDENT -- and a static term's cache is never cleared, so after the
    //! first iteration the calls are pure lookups.
    virtual void            RefreshForDensity(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override
    {
        // ONE call per term, whatever caches it owns (R1.0h + V1.35, 2026-09-10).  This used to be four
        // calls with two cross-casts and an `if constexpr`, because the real-block cache sat behind a
        // SEPARATE hook on a SEPARATE face.  Both faces now inherit `HT_BlockOwner<T>` virtually, so a term
        // with two caches presents ONE `BuildBlocks` that prepares both -- the fold stopped needing to know
        // that two caches are even possible.
        if (cd)
        {
            for (const auto& t : itsDHTs) t->RefreshForDensity(cd);           // each term's own state, for cd...
            for (const auto& t : itsHF_HTs) t->RefreshForDensity(bs, cd);     // whole-system J/K (V1.41 stage 2)
        }
        for (const auto& t : itsSHTs) t->BuildBlocks(bs, cd);
        for (const auto& t : itsDHTs) t->BuildBlocks(bs, cd);   // ...then every block from it, STAMPED: GetMatrix only reads
    }
    virtual std::ostream&   Write(std::ostream&) const override;

protected:
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::vector<std::unique_ptr<   tStatic_HT<T>>> shtv_t;
    typedef std::vector<std::unique_ptr<  tDynamic_HT<T>>> dhtv_t;
    typedef std::vector<std::unique_ptr<tDynamic_HF_HT<T>>> hf_htv_t;

    shtv_t   itsSHTs;
    dhtv_t   itsDHTs;
    const tRealSpaceResidual<T>* itsRealSpaceTerm=nullptr;   //!< the term owning the real-space quadrature, one of itsDHTs; null if none
    hf_htv_t itsHF_HTs;   // HF capable terms require a widened interface for efficient J/K table handling.

    SpinGroup itsSpinGroup;
    bool   itsIsRelativistic;
    bool   itsIsVirialValid=true;   //!< AND of the terms' IsVirialValid() (see Add); true for an empty H
    bool   itsPreservesReal=true;   //!< AND of the terms' PreservesReal() (see Add); true for an empty H
};

// NO InsertStandardTerms.  There is no "standard" term set to insert (R2.8; user 2026-08-07: the whole
// idea was too clever by half).  Every concrete Hamiltonian now spells its own core out in three Adds,
// because the core genuinely varies along axes the one helper could not express:
//   * BARE vs PSEUDISED nuclei decides BOTH the ion charge (Z vs Zion) AND the electron-nuclear term(s) --
//     and that is NOT the double/dcmplx axis it was filed under: Ham_PP is <double> and never called the
//     helper either.
//   * the BASIS lineage then decides WHICH electron-nuclear implementation (PP_Local's closed form or mesh quadrature vs
//     Ven_PP_Short/_Long's G-space route) -- a second, independent axis.
//   * the Dirac Hamiltonians replace the kinetic term outright (DiracKinetic + RestMass) and carry no
//     ion-ion at all.
// Seven of eleven Hamiltonians shared the helper; for the other four it was noise, and its <dcmplx>
// specialization was a bare assert(false).  Three explicit Adds read better than a name that hides them.
using rHamiltonianImp = tHamiltonianImp<double>;
using cHamiltonianImp = tHamiltonianImp<dcmplx>;

} //namespace