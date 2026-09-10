// File: IrrepCD.C  Exact charged density for ONE irreducable representation basis set.
//
// STEP 3c-2b (doc/RealComplexPlan.md): the leaf family is split by LINEAGE, not scalar.  V1.6/V1.7
// attached the finite faces (AO projection, HF pair partner) to T==double and the periodic face (the
// reciprocal trio) to T==dcmplx -- an encoding that quietly identified SCALAR with LINEAGE.  A real
// TRIM block's density is <double> AND periodic, which breaks that identification in both directions:
// it needs the Fourier trio and must NOT carry the AO/HF faces.  So lineage is now a CLASS axis:
//
//   IrrepCD_Core<T>      -- the shared density-matrix machinery (contractions, mixing, rho(r), ...)
//   FiniteIrrepCD        -- the FINITE leaf (molecules/atoms): core + AO projection + HF pair partner.
//   PeriodicIrrepCD<T>   -- the PERIODIC leaf (Bloch blocks): core + the reciprocal trio, and for
//                           T==double the run-typed energy-contract capability (RealBlockEnergy_CD).
//
// The scalar-keyed conditional bases are GONE, and with them the finite leaf's asserting HF stubs (the
// R2.8 denial smell): nothing instantiates a finite complex leaf, so nothing has to deny anything.
// The lineage choice is made ONCE, in IrrepCD_Factory, by probing the basis for the G-space capability.
//
// V1.32: the finite leaf is NOT a template.  It had exactly one instantiation and the factory's
// if-constexpr made a finite-complex density unrepresentable, so the parameter was vestigial -- and the
// name said the wrong thing besides.  *Finite* is this leaf's identity; the scalar is not.  (The same
// question was ASKED of PeriodicIrrepCD<T> and DECLINED: its T is load-bearing -- real TRIM block vs
// general k -- which is the whole point of the lineage-as-class split above.)
module;
#include <iosfwd>
#include <cstddef>
#include <map>
#include <string>
#include <type_traits>
export module qchem.ChargeDensity.Imp.IrrepCD;

export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;
export import qchem.ChargeDensity.FourierDensity;   // G-space rho-tilde (periodic lineage)
import qchem.ChargeDensity.Types;

export namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  This maintains the exact charge density represented by the density matrix
//  of one irreducable representation.  The full charge density will in general
//  be a summation of these guys.
//

//! \brief The reciprocal-space trio for a PERIODIC leaf block -- BOTH scalars (Step 3c-2b: a real TRIM
//! block's D contributes to ρ̃ exactly like a complex block's; the tensors follow TFit==dcmplx either
//! way, so ONE mixin body serves both).  Same CRTP shape as \c IrrepCD_HFPair: the implementation
//! reaches the block's own D and basis without any of it becoming public.
template <class Leaf> class IrrepCD_Fourier : public virtual FourierDensity
{
public:
    //! Metric-free rho-tilde(Delta-m) of THIS block: D contracted against the basis's D-free OVERLAP tensor
    //! Overlap3C(c) (empty kernel) -- the periodic density's native representation.
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // raw rho_DM (0.5(f2)); empty if no raw path
    //! V_H(Delta-m) of THIS block: D contracted against the basis's D-free Coulomb tensor Repulsion3C(c).
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;
private:
    const Leaf& self() const {return static_cast<const Leaf&>(*this);}
};

//! \brief The bra-ket PAIR-PARTNER implementation for the FINITE leaf (V1.6/V1.8).  Exact exchange is
//! finite-molecular by construction; the periodic leaf simply does not inherit this, so no denial is
//! ever written.  CRTP on the leaf: the implementation uses the block's OWN density and basis without
//! any of it becoming public (\c IrrepCD grants friendship; nothing is exposed through an accessor).
template <class Leaf> class IrrepCD_HFPair : public virtual tHF_Pair_CD<double>
{
public:
    virtual void AccumulateDirectBoth  (rsmat_t& Ji, rsmat_t& Jj, const tHF_Pair_CD<double>& other) const;
    virtual void AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const tHF_Pair_CD<double>& other) const;
    virtual void CompleteDirectPair    (rsmat_t& Ji, rsmat_t& Jj, const rsmat_t& Di,
                                        const BasisSet::Orbital_HF_IBS<double>* bs_i) const;
    virtual void CompleteExchangePair  (rsmat_t& Ki, rsmat_t& Kj, const rsmat_t& Di,
                                        const BasisSet::Orbital_HF_IBS<double>* bs_i) const;
private:
    const Leaf& self() const {return static_cast<const Leaf&>(*this);}
};

//! \brief THE SHARED DENSITY-MATRIX CORE: every operation both lineages support identically -- the
//! contractions, mixing/convergence, rho(r), the freshness serial.  Lineage-specific capabilities live
//! on the two leaves below; this class carries NO lineage assumption at all.
template <class T> class IrrepCD_Core
    : public virtual tDM_CD<T>
{
public:
    typedef  mat_t<T>  DenMat;
    typedef hmat_t<T> DenSMat; //Density matrix: HERMITIAN (= symmetric for real T, byte-identical there).
    using scalar_t = T;        //!< the block scalar (the Fourier mixin's basis cast is keyed on it)

    IrrepCD_Core();
    IrrepCD_Core(const DenSMat&,const tobs_t<T>*, Irrep);

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;
    virtual rvec_t ProjectOnto(const Fitting::ScalarProjector&) const;
    virtual double GetTotalCharge(                      ) const;

    virtual size_t Version() const {return itsVersion;}

    virtual void   ReScale      (double factor              )      ; // No UT coverage
    virtual void   MixIn        (const tMixableDensity<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;

    virtual std::ostream&       Write(std::ostream&) const;

protected:
    bool IsZero() const;

    DenSMat          itsDensityMatrix;
    const tobs_t<T>* itsBasisSet;
    Spin             itsSpin;
    Irrep            itsIrrep;
    size_t           itsVersion;   //!< TRANSIENT freshness serial (NextDensityVersion); never serialize.
};

//! \brief THE FINITE LEAF (molecules/atoms): the core plus the finite-only capabilities -- the AO
//! (auxiliary-basis) projection, the exact-exchange pair partner, and the whole-system Fock route.
//! NOT a template (V1.32): a finite density is real, full stop, so it names its bases outright rather
//! than reaching them through scalar-keyed conditionals that had one live branch each.
class FiniteIrrepCD
    : public IrrepCD_Core<double>
    , public Fitting::CoulombMetric_ProjectedDensity   // AO projection (was ProjectedDensityBase<T>)
    , public IrrepCD_HFPair<FiniteIrrepCD>             // exact-exchange pair partner (was IrrepHF_PairBase)
{
public:
    using IrrepCD_Core<double>::IrrepCD_Core;   // the core's ctors are the leaf's

    //! V1.31 whole-system route: this block's basis answers the capability, and the block folds its own
    //! density up to AO.
    virtual const BasisSet::WholeSystemFock_IBS<double>* WholeSystemFock() const;
    virtual void AddAODensity(rsmat_t& Dao) const;
    //! AO (auxiliary-basis) projection <rho|c> -- the finite path's ProjectedDensity_AO face.
    virtual double FitGetConstraint() const {return this->GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;
    //! \f$\nabla\rho\f$ from the density matrix (the molecular contraction).
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

private:
    friend class IrrepCD_HFPair<FiniteIrrepCD>;   // uses this block's own D/basis; nothing is exposed publicly
    //! The diagonal (self-paired) HF contraction; called only from the pair mixin's self-pair branch.
    void AccumulateDirect  (rsmat_t& Jii) const;
    void AccumulateExchange(rsmat_t& Kii) const;
};

//! Conditional real-block-energy base for the periodic leaf: ONLY the \c <double> instantiation (a real
//! block inside a complex run) carries the run-typed energy-contract capability.  This conditional keys
//! on exactly what it means -- the SCALAR of a periodic leaf -- unlike the retired lineage-by-scalar ones.
struct NoRealBlockEnergy {};
template <class T> using PeriodicRealEnergyBase =
    std::conditional_t<std::is_same_v<T,double>, RealBlockEnergy_CD, NoRealBlockEnergy>;

//! \brief THE PERIODIC LEAF (Bloch blocks, both scalars -- Step 3c-2b): the core plus the reciprocal
//! trio.  Carries NEITHER the AO projection NOR the HF faces (no periodic exact exchange, no auxiliary-
//! Gaussian fit), so a cross-cast probe on it tells the truth.  The \c <double> instantiation is the
//! real TRIM block's density -- real D, real DM-side GEMMs -- and additionally answers the composite's
//! run-typed energy contraction (\c RealBlockEnergy_CD).
template <class T> class PeriodicIrrepCD
    : public IrrepCD_Core<T>
    , public IrrepCD_Fourier<PeriodicIrrepCD<T>>   // the reciprocal trio (both scalars)
    , public PeriodicRealEnergyBase<T>             // run-typed energy contract (double only)
{
public:
    using IrrepCD_Core<T>::IrrepCD_Core;

    //! Periodic overlap is uncached-complex or cached-real by scalar; one body, if-constexpr split.
    virtual double GetTotalCharge() const;
    //! The periodic path is LDA-only: grad(rho) is not wired, and a silent zero would hand a GGA a
    //! plausible wrong field (R1.4) -- THROW.
    virtual rvec3_t Gradient(const rvec3_t&) const;
    //! The run-typed energy contraction (RealBlockEnergy_CD; meaningful for T==double only -- on the
    //! complex instantiation this overrides nothing and is never called).
    virtual double DM_ContractE(const Dynamic_CC_RealBlock*, const tChargeDensity<dcmplx>*) const;

private:
    friend class IrrepCD_Fourier<PeriodicIrrepCD<T>>;  // the trio uses this block's own D/basis
};

//! \brief SAME DENSITY, SAME VALUE, CHEAPER ROUTE: D stays the truth on the \a Leaf; only \f$\rho(r)\f$ is
//! factored.  \f$\rho_g=\Phi_g^\dagger D\Phi_g\f$ is a sum over PAIRS costing O(npts n^2); D is a density
//! matrix, hence PSD with rank = the occupied count, so \f$D=LL^\dagger\f$ turns it into a sum over SINGLES,
//! \f$\rho_g=\lVert L^\dagger\Phi_g\rVert^2\f$, costing O(npts n r).  Measured r=14-19 against n=118 on the
//! MnO benchmark (\c GPW_DM_RANK=1).  Exact to roundoff -- a cost change, not an accuracy trade.
//!
//! WHY A DERIVED LEAF AND NOT A DENSITY TYPE (doc/OpenWork.md, the design ruling).  Low rank is NOT CLOSED
//! UNDER MIXING: \f$(1-c)L_1L_1^\dagger+cL_2L_2^\dagger\f$ has no rank-r factor, so a first-class factored
//! density would satisfy \c tMixableDensity's SIGNATURE while breaking its behaviour -- silently, as rank
//! creep.  And of \c tDM_CD's four operations the factored form wins ONE: the two contractions go from
//! O(n^2) to O(n^2 r).  So the factor is a DERIVED, CACHED representation used only by the operation that
//! benefits; D, \c MixIn and the contractions are INHERITED UNCHANGED and LSP holds by construction.
//!
//! ONE TEMPLATE, so the two orthogonal axes (leaf x factorisation) COMPOSE instead of multiplying: the
//! factory picks the leaf from the basis and the route from its argument, and \c using \c Leaf::Leaf makes
//! every combination constructible identically.
template <class Leaf> class FactoredRho : public Leaf
{
    using T = typename Leaf::scalar_t;
public:
    using Leaf::Leaf;                       // as the leaves already do from the core

    //! \copydoc IrrepCD_Core::ProjectOnto
    //! Factored route; falls back to \c Leaf::ProjectOnto (the full quadratic form) whenever the factor
    //! does not apply -- D not PSD, or a rank too fat to pay for itself.  Same values either way.
    virtual rvec_t ProjectOnto(const Fitting::ScalarProjector&) const;

private:
    //! The factor, memoized per DENSITY SERIAL.  \c ProjectOnto is const and the factorisation is a pure
    //! function of D, so this is a normal derived cache; the invalidation key already exists, because
    //! \c IrrepCD_Core bumps \c itsVersion on every mutation of D (\c ReScale and \c MixIn are the only two).
    mutable mat_t<T> itsL;
    mutable size_t   itsRank          = 0;
    mutable size_t   itsFactorVersion = size_t(-1);   //!< the IrrepCD_Core::Version() this factor was built from
    mutable bool     itsFactorable    = false;        //!< false => the fallback is in force for this serial
    //! Tr(D) at factorisation time.  A STALE MEMO IS A SILENTLY WRONG rho -- the failure mode with no
    //! symptom -- so the version key is not trusted alone: this O(n) trace is re-checked on every memo hit.
    mutable double   itsFactorTrace   = 0.0;
};

} //namespace
