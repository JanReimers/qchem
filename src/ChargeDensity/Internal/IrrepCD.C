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
#include <memory>
#include <vector>   // JoinFold (the k-fold)
export module qchem.ChargeDensity.Imp.IrrepCD;
import qchem.Structure;   // DM_Pullback (GeometryOptimizationPlan §13.2)

export import qchem.ChargeDensity;
export import qchem.ChargeDensity.DensityMatrixDistance;
export import qchem.Symmetry.Irrep;
import qchem.ChargeDensity.Internal.FoldableBlock;   // the k-fold face a periodic block answers (not re-exported: Internal)
import qchem.ChargeDensity.Types;

export namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  This maintains the exact charge density represented by the density matrix
//  of one irreducable representation.  The full charge density will in general
//  be a summation of these guys.
//

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
    virtual rvec3vec_t DirectPullbackWith      (const tHF_Pair_CD<double>& other, const Structure* st) const;
    virtual rvec3vec_t ExchangePullbackWith    (const tHF_Pair_CD<double>& other, const Structure* st) const;
    virtual rvec3vec_t CompleteDirectPullback  (const rsmat_t& Di, const BasisSet::Orbital_HF_IBS<double>* bs_i, const Structure* st) const;
    virtual rvec3vec_t CompleteExchangePullback(const rsmat_t& Di, const BasisSet::Orbital_HF_IBS<double>* bs_i, const Structure* st) const;
private:
    const Leaf& self() const {return static_cast<const Leaf&>(*this);}
    rvec3vec_t CompletePullback(const rsmat_t& Di, const BasisSet::Orbital_HF_IBS<double>* bs_i, const Structure* st, bool exchange) const;
};

//! \brief THE SHARED DENSITY-MATRIX CORE: every operation both lineages support identically -- the
//! contractions, mixing/convergence, rho(r), the freshness serial.  Lineage-specific capabilities live
//! on the two leaves below; this class carries NO lineage assumption at all.
template <class T> class IrrepCD_Core
    : public virtual tDM_CD<T>
    , public virtual tDensityMatrixDistance<T>
    , public virtual FoldableBlock
{
public:
    typedef  mat_t<T>  DenMat;
    typedef hmat_t<T> DenSMat; //Density matrix: HERMITIAN (= symmetric for real T, byte-identical there).
    using scalar_t = T;        //!< the block scalar (the Fourier mixin's basis cast is keyed on it)

    IrrepCD_Core();
    IrrepCD_Core(DenSMat,const tobs_t<T>*, Irrep);   //!< D BY VALUE: an rvalue is moved in (a density is n×n)

    virtual double DM_Contract(const tStatic_CC<T>*) const override;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const override;
    virtual rvec3vec_t DM_Pullback(const tStatic_PC<T>* v, const Structure* st) const override {return v->Pullback(itsBasisSet, itsSpin, itsDensityMatrix, st);}
    virtual rmat3d_t   DM_Virial  (const tStatic_VC<T>* v, const Structure* st) const override {return v->Virial  (itsBasisSet, itsSpin, itsDensityMatrix, st);}
    virtual double DM_ContractBlocks(const std::map<Irrep,hmat_t<T>>&) const override;
    //! \copydoc tChargeDensity::Project -- the PERIODIC lineage: contract my D into the block's forward.
    virtual Projections Project(const BasisSet::PairProjector<dcmplx>&) const override;
    using tChargeDensity<T>::Project;   // the molecular overload (FiniteIrrepCD overrides it)
    //! \copydoc FoldableBlock::JoinFold -- a non-zero block joins with its orbital block and its D (its \c Project contracts exactly that).
    virtual bool JoinFold(const BasisSet::PairProjector<dcmplx>& p, std::vector<BasisSet::FoldBlockRef>& blocks,
                          std::vector<BasisSet::FoldBlockD>& Ds) const override;
    //! \copydoc FoldableBlock::JoinPullback -- my \c DM_Pullback / \c DM_Virial are exactly the term's per-block call on (basis, spin, D).
    virtual bool JoinPullback(DensityBlocks<double>& b) const override;
    virtual bool JoinPullback(DensityBlocks<dcmplx>& b) const override;
    //! The old face, kept one commit (V-ONE-OF-N 3a): forwards to \c Project.
    virtual double GetTotalCharge(                      ) const override;

    virtual size_t Version() const override {return itsVersion;}
    //! THE one density that IS \f$\rho[D]\f$ and nothing else (pin 28, user ruling 2026-11-15): the only true answer in the tree.
    virtual bool   HasNoFitComponent() const override {return true;}

    virtual void   ReScale      (double factor              ) override      ; // No UT coverage
    virtual void   Axpy         (double c, const tMixableDensity<T>&) override   ;  //this += c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const override;  //Frobenius(delta density matrix)
    virtual double GetMaxChangeFrom(const tMixableDensity<T>&    ) const override;  //MaxAbs(delta density matrix) -- CP2K's EPS_SCF measure
    virtual double SquaredDistanceTo(const tMixableDensity<T>&   ) const override;  //Tr(dD S dD S)
    virtual double Dot          (const tMixableDensity<T>&       ) const override;  //Re Tr(this^dagger that)

    virtual double operator()(const rvec3_t&) const override;

    virtual std::ostream&       Write(std::ostream&) const;

protected:
    bool IsZero() const;

    DenSMat          itsDensityMatrix;
    const tobs_t<T>* itsBasisSet;
    Spin             itsSpin;
    Irrep            itsIrrep;
    size_t           itsVersion;   //!< TRANSIENT freshness serial (NextDensityVersion); never serialize.
    //! The \c [DM subspace] census's memory (\c QCHEM_DIAGNOSTICS=dm_rank only; empty otherwise): this block's kept eigenvectors of D
    //! and the serial they were taken at, so the census compares one SERIAL with the previous one of the SAME block.  It was a
    //! process-wide static keyed by irrep (D-MUTABLE-STATICS): it outlived a run, and a second projection of one serial compared D
    //! with itself.  Diagnostic state, not a cache: nothing but the printed line reads it.
    mutable size_t   itsModesVersion = size_t(-1);
    mutable mat_t<T> itsModes;
};

//! \brief THE FINITE LEAF (molecules/atoms): the core plus the finite-only capabilities -- the AO
//! (auxiliary-basis) projection, the exact-exchange pair partner, and the whole-system Fock route.
//! NOT a template (V1.32): a finite density is real, full stop, so it names its bases outright rather
//! than reaching them through scalar-keyed conditionals that had one live branch each.
class FiniteIrrepCD
    : public IrrepCD_Core<double>
    , public IrrepCD_HFPair<FiniteIrrepCD>             // exact-exchange pair partner (was IrrepHF_PairBase)
{
public:
    using IrrepCD_Core<double>::IrrepCD_Core;   // the core's ctors are the leaf's
    virtual std::unique_ptr<tMixableDensity<double>> MakeCopy() const override
        {return std::make_unique<FiniteIrrepCD>(this->itsDensityMatrix, this->itsBasisSet, this->itsIrrep);}

    //! V1.31 whole-system route: this block's basis answers the capability, and the block folds its own
    //! density up to AO.
    virtual const BasisSet::WholeSystemFock_IBS<double>* WholeSystemFock() const override;
    virtual void AddAODensity(rsmat_t& Dao) const override;
    //! \copydoc tChargeDensity::Project -- the MOLECULAR lineage: contract my D into the block's forward.
    virtual Projections Project(const BasisSet::PairProjector<double>&) const override;
    using IrrepCD_Core<double>::Project;   // the periodic overload
    //! \f$\nabla\rho\f$ from the density matrix (the molecular contraction).
    virtual rvec3_t  Gradient  (const rvec3_t&) const override; // No UT coverage

private:
    friend class IrrepCD_HFPair<FiniteIrrepCD>;   // uses this block's own D/basis; nothing is exposed publicly
    //! The diagonal (self-paired) HF contraction; called only from the pair mixin's self-pair branch.
    void AccumulateDirect  (rsmat_t& Jii) const;
    void AccumulateExchange(rsmat_t& Kii) const;
};

//! \brief The run-typed energy contraction (\c RealBlockEnergy_CD) of a REAL periodic block inside a complex run: this
//! block's real D against the complex run's dynamic term, through the term's real block face for the run's density.
//! CRTP on the leaf (\a Leaf, whose D and basis it reads), the shape of \c Composite_HFSystem (CompositeCD.C): only the
//! \c <double> leaf carries the capability, so only it carries the implementation, and the complex leaf declares
//! nothing it cannot override (it once declared this member unconditionally, overriding nothing and throwing on <dcmplx>).
//! The conditional keys on exactly what it means -- the SCALAR of a periodic leaf -- unlike the retired lineage-by-scalar ones.
template <class Leaf> class PeriodicIrrepCD_RealEnergy : public RealBlockEnergy_CD
{
public:
    virtual double DM_ContractE(const tDynamic_CC<double>*, DensityID) const override;
private:
    const Leaf& self() const {return static_cast<const Leaf&>(*this);}
};
struct NoRealBlockEnergy {};
template <class T, class Leaf> using PeriodicRealEnergyBase =
    std::conditional_t<std::is_same_v<T,double>, PeriodicIrrepCD_RealEnergy<Leaf>, NoRealBlockEnergy>;

//! \brief THE PERIODIC LEAF (Bloch blocks, both scalars -- Step 3c-2b): the core plus the reciprocal
//! trio.  Carries NEITHER the AO projection NOR the HF faces (no periodic exact exchange, no auxiliary-
//! Gaussian fit), so a cross-cast probe on it tells the truth.  The \c <double> instantiation is the
//! real TRIM block's density -- real D, real DM-side GEMMs -- and additionally answers the composite's
//! run-typed energy contraction (\c RealBlockEnergy_CD).
template <class T> class PeriodicIrrepCD
    : public IrrepCD_Core<T>
    , public PeriodicRealEnergyBase<T,PeriodicIrrepCD<T>>   // run-typed energy contract (double only)
{
public:
    using IrrepCD_Core<T>::IrrepCD_Core;
    //! A plain periodic leaf, also for a \c FactoredRho one: a history copy is only ever Axpy'd / Dot'd, never sampled.
    virtual std::unique_ptr<tMixableDensity<T>> MakeCopy() const
        {return std::make_unique<PeriodicIrrepCD<T>>(this->itsDensityMatrix, this->itsBasisSet, this->itsIrrep);}

    //! Periodic overlap is uncached-complex or cached-real by scalar; one body, if-constexpr split.
    virtual double GetTotalCharge() const;
    //! The periodic path is LDA-only: grad(rho) is not wired, and a silent zero would hand a GGA a
    //! plausible wrong field (R1.4) -- THROW.
    virtual rvec3_t Gradient(const rvec3_t&) const;
    // The run-typed energy contraction lives in PeriodicIrrepCD_RealEnergy (a base of the <double> leaf only).

private:
    friend class PeriodicIrrepCD_RealEnergy<PeriodicIrrepCD<T>>;  // so does the real-block energy contraction
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

    //! A block that would take its FACTORED path declines the fold (its projection is its own call); otherwise the leaf's rule.
    virtual bool JoinFold(const BasisSet::PairProjector<dcmplx>& p, std::vector<BasisSet::FoldBlockRef>& blocks,
                          std::vector<BasisSet::FoldBlockD>& Ds) const override;
    //! \copydoc IrrepCD_Core::Project
    //! Factored route; falls back to \c Leaf::Project (the full quadratic form) whenever the factor
    //! does not apply -- D not PSD, or a rank too fat to pay for itself.  Same values either way.
    virtual Projections Project(const BasisSet::PairProjector<dcmplx>&) const override;
    using Leaf::Project;

private:
    //! The factor, memoized per DENSITY SERIAL.  \c ProjectOnto is const and the factorisation is a pure
    //! function of D, so this is a normal derived cache; the invalidation key already exists, because
    //! \c IrrepCD_Core bumps \c itsVersion on every mutation of D (\c ReScale and \c MixIn are the only two).
    mutable mat_t<T> itsL;
    mutable size_t   itsRank          = 0;
    mutable size_t   itsFactorVersion = size_t(-1);   //!< the IrrepCD_Core::Version() this factor was built from
    mutable bool     itsFactorable    = false;        //!< false => the fallback is in force for this serial
    //! true while this block's D stays outside the PSD cone after its warning was printed: the canary prints once when a block LEAVES
    //! the cone and again only after it has come back (D-MUTABLE-STATICS: it was a process-wide cap of 5 warnings, silent for every
    //! later run in the process).
    mutable bool     itsNotPSDWarned  = false;
    //! Tr(D) at factorisation time.  A STALE MEMO IS A SILENTLY WRONG rho -- the failure mode with no
    //! symptom -- so the version key is not trusted alone: this O(n) trace is re-checked on every memo hit.
    mutable double   itsFactorTrace   = 0.0;
};

} //namespace
