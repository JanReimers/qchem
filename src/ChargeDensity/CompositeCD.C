// File: CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <vector>
#include <memory>
#include <cstddef>
#include <map>
#include <string>
#include <variant>
export module qchem.CompositeCD;
export import qchem.ChargeDensity;
export import qchem.ChargeDensity.FourierDensity;   // G-space rho-tilde (summed over k-blocks)
export import qchem.Symmetry.Lattice_3D.SpaceGroup; // ReciprocalOp {U|τ} -- the IBZ density-symmetrization ops (glide phase)
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{

//! \brief THE CHILD SLOT (doc/RealComplexPlan.md §4, Step 2): one block density, typed by ITS OWN scalar
//! rather than the composite's face -- so children of one composite may differ (a real TRIM block beside
//! general-k complex blocks on a mixed mesh).  The alternatives are the ABSTRACT per-block face, ownership
//! included.  Aggregation stays single-source: scalar-independent operations visit with ONE generic lambda;
//! the T-typed operations forward to the same-scalar alternative (the cross arm becomes reachable -- and
//! gets its narrowing implementation -- when Step 3 un-pins the basis type per block).
using cd_child_t = std::variant<std::unique_ptr<tDM_CD<double>>, std::unique_ptr<tDM_CD<dcmplx>>>;

//--------------------------------------------------------------------------
//
//  Full charge density represented Compositely as sum of density matricies.
//  Templated on the matrix element type T (rX/cX); the <double> alias preserves the existing
//  real callers, the <dcmplx> instantiation aggregates the plane-wave (Bloch-irrep) densities.
//
//! \brief The reciprocal-space trio for a COMPOSITE -- periodic path only (V1.7).  Each block already
//! carries its BZ weight, so every one of these is a straight sum over the leaves.  CRTP for the same reason
//! as the leaf's: reach the contained blocks without exposing them.
template <class Comp> class Composite_Fourier : public virtual FourierDensity
{
public:
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // empty if any block lacks the raw path
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;
private:
    const Comp& self() const {return static_cast<const Comp&>(*this);}
};

template <class T, class Comp> using CompositeFourierBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Composite_Fourier<Comp>, NoFourierDensity>;

//! \brief The whole-system exact-exchange sweep for a COMPOSITE -- real path only (V1.6, completing it).
//! Same CRTP shape as \c Composite_Fourier, and for the same reason the LEAF got one: inheriting the
//! \c tHF_System_CD face conditionally but DECLARING its two methods unconditionally left
//! \c tComposite_CD<dcmplx> carrying exact-exchange members that override nothing and can only throw --
//! the very "denial the other half of the hierarchy must write" this item set out to remove, surviving
//! one level above the leaf.  Declared here, the complex composite grows nothing at all.
template <class Comp> class Composite_HFSystem : public virtual tHF_System_CD<double>
{
public:
    virtual void AccumulateDirectAll  (std::vector<hmat_t<double>>& Jall) const;
    virtual void AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const;
private:
    const Comp& self() const {return static_cast<const Comp&>(*this);}
};

template <class T, class Comp> using CompositeHFBase =
    std::conditional_t<std::is_same_v<T,double>, Composite_HFSystem<Comp>, NoHF_System>;

template <class T> class tComposite_CD
    : public virtual tDM_CD<T>
    , public virtual tLineageTracked<T> // Layer-2: this top-level density tracks its SCF lineage head
    , public ProjectedDensityBase<T> // AO projection on the finite (double) path; empty on the periodic path
    , public CompositeFourierBase<T,tComposite_CD<T>>   // reciprocal trio: periodic path only (V1.7)
    , public CompositeHFBase<T,tComposite_CD<T>>        // whole-system exact exchange: real path only (V1.6)
{
public:
    //! \a pointOps = the reciprocal point group for IBZ density symmetrization (doc/GPWPlan1.md item 3).  The
    //! G-space density accessors then return the STAR AVERAGE, which is what makes an IBZ-reduced density exact
    //! (the star weights alone give only the correct band sum).  Default {} = trivial group {E} = exact no-op
    //! -- molecules / Γ / unreduced crystals pass through untouched (the general form; "no symmetry" = trivial).
    //! It is a ctor argument, not a setter: the symmetry is a fixed property of the density, set once at build.
    explicit tComposite_CD(std::vector<Symmetry::Lattice_3D::ReciprocalOp> pointOps = {});
    //! TAKES OWNERSHIP of \a cd (V1.25).  TWO overloads, one per child scalar (doc/RealComplexPlan.md
    //! Step 2): the child slot is typed by the BLOCK, not by this composite's face, so either alternative
    //! may be inserted regardless of T.  A molecular composite simply never receives the complex one.
    void Insert(std::unique_ptr<tDM_CD<double>> cd);
    void Insert(std::unique_ptr<tDM_CD<dcmplx>> cd);

    // The whole-system J/K sweep is NOT declared here (V1.6 ISP): it lives in Composite_HFSystem, which
    // only the REAL instantiation inherits -- so tComposite_CD<dcmplx> declares nothing and defines nothing
    // for an operation the periodic path never has (Ham_PW_DFT adds Vee_Hartree, never exact exchange).
    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;   // sum over irrep blocks
    virtual rvec_t DM_RhoAtPoints(const BasisSet::cFIT_SF_Delta&) const;   // sum over irrep blocks
    //! The mixed-run overload (3c-3): a REAL child GEMMs its own PhiR table instead of the pointwise fallback.

    virtual double GetTotalCharge      (                     ) const;

    // The blocks are mutated together (MixIn/ReScale fan out to all), so any block's serial tracks the
    // composite's freshness; forward to the first.  Empty composite -> 0 (the "no density yet" sentinel).
    virtual size_t Version() const
    {return itsCDs.empty() ? 0 : std::visit([](const auto& c){return c->Version();}, itsCDs.front());}

    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;

    virtual void   ReScale      (double factor         )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const tMixableDensity<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const;

    // The reciprocal trio is NOT declared here (V1.7 ISP) -- see Composite_Fourier below, inherited only
    // on the periodic path, so the finite composite is not asked questions it cannot answer.

private:
    friend class Composite_Fourier<tComposite_CD<T>>;   // sums the blocks; nothing is exposed publicly
    friend class Composite_HFSystem<tComposite_CD<T>>;  // drives the canonical-pair sweep over the blocks
    tComposite_CD(const tComposite_CD&);

    typedef std::vector<cd_child_t> cdv_t;   // the §4 child slot: per-block scalar, face-T-independent
    cdv_t itsCDs;
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> itsPointOps;   //!< reciprocal {U|τ} point group for IBZ symmetrization ({E} when empty)
};

using rComposite_CD = tComposite_CD<double>;   using cComposite_CD = tComposite_CD<dcmplx>;

} //namespace