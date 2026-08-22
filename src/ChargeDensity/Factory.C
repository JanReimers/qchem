// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
export module qchem.ChargeDensity.Factory;
export import qchem.ChargeDensity;
export import qchem.FittedCD;
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{
    typedef std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs_t;   //!< the Coulomb-metric (density-fit) face

    //! \brief WHICH ROUTE evaluates \f$\rho(r)\f$ from D -- the ONLY factorisation detail that crosses this
    //! boundary (doc/OpenWork.md, the factored-rho design ruling).  D stays the truth in every case: the
    //! routes differ in COST, not in value (exact to roundoff), and a route that cannot apply falls back to
    //! \c Direct by INHERITANCE rather than by a branch.
    //!
    //! The LEAF axis is deliberately NOT an enum: the factory derives it from the basis it is handed, so a
    //! caller cannot request a leaf inconsistent with its own argument.  One enum, two orthogonal axes.
    enum class RhoRoute
    {
        Direct,           //!< \f$\rho_g=\Phi_g^\dagger D\,\Phi_g\f$ -- O(npts n^2); the leaf's own body.
        PivotedCholesky,  //!< \f$D=LL^\dagger\f$, \f$\rho_g=\lVert L^\dagger\Phi_g\rVert^2\f$ -- O(npts n r).
        EigenTrim,        //!< NOT WIRED -- throws.  The minimal-rank factor; the route for a D that has left
                          //!< the PSD cone, where Cholesky is inapplicable.  Listed so the gap is visible.
    };
    //! The route used when a caller does not name one.  \c PivotedCholesky, unless \c QCHEM_DM_LOWRANK=0
    //! forces \c Direct -- the A/B valve, and the escape hatch if a future density ever needs it.
    RhoRoute DefaultRhoRoute();

    //! \a route defaults (via \c DefaultRhoRoute) to the factored \f$\rho\f$; \c Direct is the plain leaf.
    template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& DM,const tobs_t<T>* bs, Irrep,
                                                  RhoRoute route);
    //! Same, taking the default route -- so every existing call site is untouched.
    template <class T> tDM_CD<T>* IrrepCD_Factory(const hmat_t<T>& DM,const tobs_t<T>* bs, Irrep); // DM Hermitian
    //! Build a polarized density from its two channels, TAKING OWNERSHIP of both (V1.25).
    template <class T> std::unique_ptr<tDM_CD<T>>
        PolarizedCD_Factory(std::unique_ptr<tDM_CD<T>> up, std::unique_ptr<tDM_CD<T>> down);
    std::unique_ptr<FittedCD> FittedCD_Factory(fbs_t&, double totalCharge); //!< caller owns the result

} //namespace