// File: IrrepCD.C  Exact charged density for ONE irreducable representation basis set.
//
// STEP 3c-2b (doc/RealComplexPlan.md): the leaf family is split by LINEAGE, not scalar.  V1.6/V1.7
// attached the finite faces (AO projection, HF pair partner) to T==double and the periodic face (the
// reciprocal trio) to T==dcmplx -- an encoding that quietly identified SCALAR with LINEAGE.  A real
// TRIM block's density is <double> AND periodic, which breaks that identification in both directions:
// it needs the Fourier trio and must NOT carry the AO/HF faces.  So lineage is now a CLASS axis:
//
//   IrrepCD_Core<T>      -- the shared density-matrix machinery (contractions, mixing, rho(r), ...)
//   IrrepCD<T>           -- the FINITE leaf (molecules/atoms): core + AO projection + HF pair partner.
//                           The name is unchanged, so every molecular consumer is untouched.
//   PeriodicIrrepCD<T>   -- the PERIODIC leaf (Bloch blocks): core + the reciprocal trio, and for
//                           T==double the run-typed energy-contract capability (RealBlockEnergy_CD).
//
// The scalar-keyed conditional bases are GONE, and with them IrrepCD<dcmplx>'s asserting HF stubs (the
// R2.8 denial smell): nothing instantiates a finite complex leaf, so nothing has to deny anything.
// The lineage choice is made ONCE, in IrrepCD_Factory, by probing the basis for the G-space capability.
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
    virtual rvec_t DM_RhoAtPoints(const rvec3vec_t&, const std::map<Irrep,mat_t<T>>&) const;
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

//! The pair face on the finite leaf's double instantiation (the only one that exists), nothing otherwise.
template <class T, class Leaf> using IrrepHF_PairBase =
    std::conditional_t<std::is_same_v<T,double>, IrrepCD_HFPair<Leaf>, NoHF_Pair>;

//! \brief THE FINITE LEAF (molecules/atoms; the historical name, so molecular consumers are untouched):
//! the core plus the finite-only capabilities -- the AO (auxiliary-basis) projection, the exact-exchange
//! pair partner, and the whole-system Fock route.  Only the \c <double> instantiation exists (a finite
//! complex density is not a thing); the templated form is kept for symmetry with the family.
template <class T> class IrrepCD
    : public IrrepCD_Core<T>
    , public ProjectedDensityBase<T> // AO projection (CoulombMetric_ProjectedDensity for double)
    , public IrrepHF_PairBase<T,IrrepCD<T>>   // exact-exchange pair partner (double)
{
public:
    using IrrepCD_Core<T>::IrrepCD_Core;   // the core's ctors are the leaf's

    //! V1.31 whole-system route: this block's basis answers the capability, and the block folds its own
    //! density up to AO.
    virtual const BasisSet::WholeSystemFock_IBS<T>* WholeSystemFock() const;
    virtual void AddAODensity(hmat_t<T>& Dao) const;
    //! AO (auxiliary-basis) projection <rho|c> -- the finite path's ProjectedDensity_AO face.
    virtual double FitGetConstraint() const {return this->GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;
    //! \f$\nabla\rho\f$ from the density matrix (the molecular contraction).
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

private:
    friend class IrrepCD_HFPair<IrrepCD<T>>;   // uses this block's own D/basis; nothing is exposed publicly
    //! The diagonal (self-paired) HF contraction; called only from the pair mixin's self-pair branch.
    void AccumulateDirect  (hmat_t<T>& Jii) const;
    void AccumulateExchange(hmat_t<T>& Kii) const;
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

} //namespace
