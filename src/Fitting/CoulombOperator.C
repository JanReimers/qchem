// File: Fitting/CoulombOperator.C  The DENSITY fits a Hartree term holds: the first-order face and the Dunlap face, and their factories.
//
// (The module keeps its 2026-11-06 name; its operators moved to qcMesh in V-ONE-OF-N 3b-ii, and since 3c it is the home of the
// density-fit faces the qchem.Fitting.FunctionFitter hub used to carry.)  A fit basis says what a function is and answers its own
// Coulomb integrals (ruling (f)); a fitter projects, solves in its own metric, and keeps its coefficients private (ruling (b)).
// The concrete classes stay Internal: a caller holds these faces and the Make functions below, where the concrete pairing lives.
//
// MEASURED (lane S prototype, TeF₆ at a fixed density, doc/Records/ on the row): at a mesh's own rung the J error of any
// route that forms v at points is floored by the final contraction, 4.7e-4 per element at 80/29.  Analytic (ij|a) has no
// such floor.  The direct sum sits well above it (about 1e-2, flat in the radial count).
module;
#include <functional>
#include <memory>
export module qchem.Fitting.CoulombOperator;
export import qchem.BasisSet.PairProjector;    // what both density-fit faces are (the density projects itself onto them)
export import qchem.BasisSet.Orbital_DFT_IBS;  // the fit-basis faces and the orbital block
export import qchem.Mesh;                     // qcMesh::Mesh (the points a direct sum runs over)
export import qchem.Mesh.SiteMultipoles;      // what a radial solve runs over
export import qchem.Mesh.CoulombOperator;     // the point set's Coulomb integrals (moved to qcMesh, V-ONE-OF-N 3b-ii)

export namespace qchem::Fitting
{

//! \brief THE FIRST-ORDER FITTED HARTREE FACE, GENERIC OVER THE FIT BASIS (V-ONE-OF-N 3b; user rulings 2026-11-08 (b), (f),
//! (g), doc/Records/VOneOfN_D2Design.md §1a): what \c Vee_Fit holds, whichever fit basis -- Gaussian, \f$\delta\f$ mesh,
//! plane wave -- the fitter is over.
//!   - The density's projections onto the fit basis arrive through \c DoFit, from \c tChargeDensity::Project(*this): the fitter
//!     IS the projector.  It keeps its coefficients PRIVATE, in its own type, and does its own solve (identity on an
//!     orthogonal basis).
//!   - \c Repulsion is the PHYSICS round: the private coefficients contracted with the fit basis's COULOMB 3-centre integrals,
//!     \f$\sum_f c_f\langle ab|1/r|f\rangle\f$.  So the term applies the repulsion by asking for it, and the basis answers in its own
//!     realization (a \f$\delta\f$ mesh through its Coulomb operator; a plane wave as \f$k(G)\langle ab|G\rangle\f$).  No basis-kind
//!     branch anywhere above.
//!   - \c HartreeEnergy is the first-order energy; \a trDV is asked for only by a fitter whose energy needs it.
//! One specialization per fit lineage, as \c PairProjector: the periodic one serves both block scalars.
template <class T> class FirstOrderDensityFit;
template <> class FirstOrderDensityFit<double> : public virtual BasisSet::PairProjector<double>
{
public:
    virtual ~FirstOrderDensityFit() = default;
    virtual void DoFit(const Projections& rho)=0;                                          //!< \f$\langle f|\rho\rangle\f$ from \c Project(*this)
    virtual hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>&) const=0; //!< \f$\sum_f c_f\langle ab|1/r|f\rangle\f$
    virtual double HartreeEnergy(const std::function<double()>& trDV) const=0;             //!< the first-order \f$E_H\f$ of the last fit
    virtual double Integral() const=0;                                                     //!< the fit's charge \f$\int\tilde\rho\f$
    //! \name THE DERIVATIVE of \c Repulsion beside it, at the FIXED fit (GeometryOptimizationPlan §13.2, lane B 2026-11-09; the user's
    //! 2026-11-08 pre-approval): \f$F_A=\partial_{R_A}\mathrm{Re}\,\mathrm{Tr}(D\,\mathrm{Repulsion}(orb))\f$, the block's centres moving.
    //! On an orthonormal fit basis whose functions do not move that is the first-order Hartree term's whole fixed-D force: \f$E_H\f$ is
    //! quadratic in D through the projection, and its R-derivative at fixed D is this, once.  \c Can is false exactly when \c Make refuses.
    //! CLIENT: Vee_Fit's force.
    //!@{
    virtual bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>&) const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>&, const hmat_t<double>& D,
                                                                     const Structure* st) const=0;
    //!@}
};
template <> class FirstOrderDensityFit<dcmplx> : public virtual BasisSet::PairProjector<dcmplx>
{
public:
    virtual ~FirstOrderDensityFit() = default;
    virtual void DoFit(const Projections& rho)=0;
    virtual hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&) const=0;   //!< a real TRIM block
    virtual hmat_t<dcmplx> Repulsion(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&) const=0;   //!< a Bloch block
    virtual double HartreeEnergy(const std::function<double()>& trDV) const=0;
    virtual double Integral() const=0;
    //! \name The derivative of \c Repulsion beside it, as on \c FirstOrderDensityFit<double>: a real TRIM block and a Bloch block.
    //!@{
    virtual bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&) const=0;
    virtual bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&) const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&, const hmat_t<double>& D,
                                                                     const Structure* st) const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&, const hmat_t<dcmplx>& D,
                                                                     const Structure* st) const=0;
    //!@}
};

//! \brief THE DUNLAP FITTED HARTREE FACE, GENERIC OVER THE FIT BASIS (V-ONE-OF-N 3b, user ruling 2026-11-08 (g)): what
//! \c Vee_FitDunlap holds.  The second-order formula \f$2E_{fit}-E_{fitfit}\f$ is the TERM's; the fitter supplies its pieces.
//!   - \c DoFit takes the density's projections on the fitter's own axis (from \c tChargeDensity::Project(*this)) and the charge the
//!     fit must carry, \a charge: Dunlap's fit is charge-CONSTRAINED, and the charge is the density's, not the projections'.
//!   - \c Repulsion is the physics round, as on \c FirstOrderDensityFit; \c SelfRepulsion is
//!     \f$\tfrac12\langle\tilde\rho|1/r_{12}|\tilde\rho\rangle\f$, the fit basis's Coulomb metric.
//! Realized today on a Gaussian fit basis only (\c ConstrainedFF).  A plane-wave one is ACCEPTED by the ruling (its metric is
//! diagonal, so it reduces to \c Vee_Fit) and not built; a \f$\delta\f$ basis has no Coulomb metric to offer (OpenWork row).
template <class T> class DunlapDensityFit;
template <> class DunlapDensityFit<double> : public virtual BasisSet::PairProjector<double>
{
public:
    virtual ~DunlapDensityFit() = default;
    virtual void DoFit(const Projections& rho, double charge)=0;                           //!< fit, constrained to \a charge
    virtual hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>&) const=0; //!< \f$\sum_f c_f\langle ab|1/r|f\rangle\f$
    virtual double SelfRepulsion() const=0;                                                //!< \f$\tfrac12\langle\tilde\rho|1/r_{12}|\tilde\rho\rangle\f$
    virtual double Integral() const=0;                                                     //!< the fit's charge
    //! \name The geometry derivatives AT FIXED FIT COEFFICIENTS (GeometryOptimizationPlan §13.2, lane B; the user's 2026-11-08
    //! pre-approval of derivative faces beside a face's methods).  What the Dunlap force needs, since the constrained fit is
    //! variational and N is fixed.  CLIENT: Vee_FitDunlap's force.  Each \c Can is false exactly when its \c Make would refuse.
    //!@{
    //! \f$F_A=\sum_{ij}D_{ij}\,\partial\,\mathrm{Repulsion}(orb)_{ij}/\partial R_A\f$, the block's and the fit's centres moving.
    virtual bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>&) const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>&, const hmat_t<double>& D,
                                                                     const Structure* st) const=0;
    //! \f$\partial\,\mathrm{SelfRepulsion}()/\partial R_A\f$, the fit's centres moving -- scaled as the value is (the ½ included).
    virtual bool CanPullbackSelfRepulsion() const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeSelfRepulsionPullback(const Structure* st) const=0;
    //!@}
};

//! \brief The Dunlap fit on a Gaussian (non-orthogonal, Coulomb-metric) fit basis \a fbs.  THROWS when \a fbs has no Coulomb metric.
std::unique_ptr<DunlapDensityFit<double>> MakeDunlapDensityFitter(std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs);

//! \brief A density fit on \a deltaBasis (an orthogonal, \f$\delta\f$-like fit basis that answers its own Coulomb integrals,
//! \c BasisSet::PointCoulombIntegrals: a \c CoulombDeltaFit_IBS or a \c MultipoleFit_IBS).  THROWS when the basis is not orthogonal or
//! carries no Coulomb integrals.
std::unique_ptr<FirstOrderDensityFit<double>>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis);

//! \brief A density fit on a PLANE-WAVE fit basis \a pw (orthonormal {G}): the fit is \f$\tilde\rho(G)\f$ itself, held privately;
//! the physics is the basis's Coulomb integrals.  THROWS when \a pw is not orthonormal or lacks the G-space faces.
std::unique_ptr<FirstOrderDensityFit<dcmplx>> MakeGSpaceDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS> pw);

} //namespace
