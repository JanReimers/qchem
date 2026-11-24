// File: ChargeDensity/Factory.C  Create some charge densitytypes.
module;
#include <memory>
export module qchem.ChargeDensity.Factory;
export import qchem.ChargeDensity;
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{

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
    //! A block density from \a DM (Hermitian; BY VALUE: pass an rvalue and it is moved in, not copied) on \a route: \c PivotedCholesky is
    //! the factored \f$\rho\f$, \c Direct the plain leaf.  Every caller STATES its route; a run's blocks are made by its \c DensityMaker.
    template <class T> tDM_CD<T>* IrrepCD_Factory(hmat_t<T> DM,const tobs_t<T>* bs, Irrep,
                                                  RhoRoute route);

    //! \brief THE RUN'S DENSITY MAKER (D-MUTABLE-STATICS-RUNPOLICY group E, the user's ruling 2026-11-23, E2'): every block density a run
    //! builds from a D -- each iteration's (the orbitals'), the seed's, a restart's -- is made here, on the ONE \c RhoRoute the run's SETUP
    //! resolved (a stated policy.dmLowRank, cp2kCompat, or the qchem default) and handed down at construction.
    //!
    //! It DECIDES NOTHING (the user: "the decision could be pre determined and passed into the Density maker at construction time"): it
    //! never switches route, and the policy that chooses one lives in the setup, which reports it.  A concrete value, not an interface --
    //! one policy, one enum (the user: an interface per policy "is overkill if there just a single policy defined by one enum RhoRoute").
    //! It replaced \c DefaultRhoRoute(), which read the process-wide run policy, so a molecule after a cp2kCompat solid took the solid's
    //! route (OpenWork §3 D-RUNPOLICY-LEAK).
    class DensityMaker
    {
    public:
        explicit DensityMaker(RhoRoute route) : itsRoute(route) {}
        //! The block density of \a DM on \a bs, irrep \a qns, on this run's route.
        template <class T> std::unique_ptr<tDM_CD<T>> Make(hmat_t<T> DM, const tobs_t<T>* bs, Irrep qns) const
        {
            return std::unique_ptr<tDM_CD<T>>(IrrepCD_Factory<T>(std::move(DM), bs, std::move(qns), itsRoute));
        }
        bool operator==(const DensityMaker&) const = default;
    private:
        RhoRoute itsRoute;
    };

} //namespace