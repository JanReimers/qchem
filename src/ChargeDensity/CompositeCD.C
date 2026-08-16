// File: CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <vector>
#include <memory>
#include <cstddef>
#include <map>
#include <string>
export module qchem.CompositeCD;
export import qchem.ChargeDensity;
export import qchem.ChargeDensity.FourierDensity;   // G-space rho-tilde (summed over k-blocks)
export import qchem.Symmetry.Lattice_3D.SpaceGroup; // ReciprocalOp {U|τ} -- the IBZ density-symmetrization ops (glide phase)
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{

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

template <class T> class tComposite_CD
    : public virtual tDM_CD<T>
    , public virtual tLineageTracked<T> // Layer-2: this top-level density tracks its SCF lineage head
    , public ProjectedDensityBase<T> // AO projection on the finite (double) path; empty on the periodic path
    , public CompositeFourierBase<T,tComposite_CD<T>>   // reciprocal trio: periodic path only (V1.7)
    , public HF_SystemBase<T>        // exact exchange spans every block: real path only, empty for dcmplx (V1.6)
{
public:
    //! \a pointOps = the reciprocal point group for IBZ density symmetrization (doc/GPWPlan1.md item 3).  The
    //! G-space density accessors then return the STAR AVERAGE, which is what makes an IBZ-reduced density exact
    //! (the star weights alone give only the correct band sum).  Default {} = trivial group {E} = exact no-op
    //! -- molecules / Γ / unreduced crystals pass through untouched (the general form; "no symmetry" = trivial).
    //! It is a ctor argument, not a setter: the symmetry is a fixed property of the density, set once at build.
    explicit tComposite_CD(std::vector<Symmetry::Lattice_3D::ReciprocalOp> pointOps = {});
    //! TAKES OWNERSHIP of \a cd.  The signature used to say `tDM_CD<T>*` while the body wrapped it in a
    //! unique_ptr -- an ownership transfer visible only by reading the implementation (V1.25).
    void Insert(std::unique_ptr<tDM_CD<T>> cd);

    virtual void AccumulateDirectAll  (std::vector<hmat_t<T>>& Jall) const;
    virtual void AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall) const;

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;   // sum over irrep blocks
    virtual rvec_t DM_RhoAtPoints(const rvec3vec_t&, const std::map<Irrep,mat_t<T>>&) const;   // sum over irrep blocks

    virtual double GetTotalCharge      (                     ) const;

    // The blocks are mutated together (MixIn/ReScale fan out to all), so any block's serial tracks the
    // composite's freshness; forward to the first.  Empty composite -> 0 (the "no density yet" sentinel).
    virtual size_t Version() const {return itsCDs.empty() ? 0 : itsCDs.front()->Version();}

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
    tComposite_CD(const tComposite_CD&);

    typedef std::vector<std::unique_ptr<tDM_CD<T>>> cdv_t;
    cdv_t itsCDs;
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> itsPointOps;   //!< reciprocal {U|τ} point group for IBZ symmetrization ({E} when empty)
};

using rComposite_CD = tComposite_CD<double>;   using cComposite_CD = tComposite_CD<dcmplx>;

} //namespace