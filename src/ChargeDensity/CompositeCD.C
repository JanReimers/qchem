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
template <class T> class tComposite_CD
    : public virtual tDM_CD<T>
    , public virtual tLineageTracked<T> // Layer-2: this top-level density tracks its SCF lineage head
    , public ProjectedDensityBase<T> // AO projection on the finite (double) path; empty on the periodic path
    , public FourierDensityBase<T>   // FourierDensity on the periodic (dcmplx) path; empty on the finite path
{
public:
    //! \a pointOps = the reciprocal point group for IBZ density symmetrization (doc/GPWPlan1.md item 3).  The
    //! G-space density accessors then return the STAR AVERAGE, which is what makes an IBZ-reduced density exact
    //! (the star weights alone give only the correct band sum).  Default {} = trivial group {E} = exact no-op
    //! -- molecules / Γ / unreduced crystals pass through untouched (the general form; "no symmetry" = trivial).
    //! It is a ctor argument, not a setter: the symmetry is a fixed property of the density, set once at build.
    explicit tComposite_CD(std::vector<Symmetry::Lattice_3D::ReciprocalOp> pointOps = {});
    void Insert(tDM_CD<T>*);

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

    //! Sum the contained blocks' G-space coefficients: \f$\tilde\rho(\Delta m)=\sum_k w_k\tilde\rho_k\f$
    //! (each block already carries its BZ weight).  Plane-wave (dcmplx) path; NA-asserts for double.
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // raw rho_DM (0.5(f2)); empty if any block lacks it
    //! Sum the contained blocks' Coulomb projections \f$V_H=\sum_k w_k V_{H,k}\f$ (\f$V_H\f$ is linear in
    //! \f$\tilde\rho\f$).  Plane-wave (dcmplx) path; NA-asserts for double.
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;

private:
    tComposite_CD(const tComposite_CD&);

    typedef std::vector<std::unique_ptr<tDM_CD<T>>> cdv_t;
    cdv_t itsCDs;
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> itsPointOps;   //!< reciprocal {U|τ} point group for IBZ symmetrization ({E} when empty)
};

using rComposite_CD = tComposite_CD<double>;   using cComposite_CD = tComposite_CD<dcmplx>;

} //namespace