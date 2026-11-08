// File: Fitting/CoulombOperator.C  How a grid fit's potential is found, kept apart from what the fit IS.
//
// A fit basis says what a function is; a CoulombOperator says how the potential of a density on it is found (user ruling
// 2026-11-06, on lane S's design note).  The Gaussian route needs neither name: (ij|a) never forms v at points, so
// ConstrainedFF stays as it is.  The grid routes (a density sampled on a mesh) all end the same way -- v at the mesh's
// points, then the same quadrature XC uses -- and differ only in how v is found, so that is the one thing that varies.
//
// Its own module, beside the qchem.Fitting.FunctionFitter hub rather than on it.  The concrete classes stay Internal: a
// caller holds this face and the Make functions below, which is where the concrete pairing lives.
//
// MEASURED (lane S prototype, TeF₆ at a fixed density, doc/Records/ on the row): at a mesh's own rung the J error of any
// route that forms v at points is floored by the final contraction, 4.7e-4 per element at 80/29.  Analytic (ij|a) has no
// such floor.  The direct sum sits well above it (about 1e-2, flat in the radial count).
module;
#include <functional>
#include <memory>
export module qchem.Fitting.CoulombOperator;
export import qchem.Fitting.FunctionFitter;   // FunctionFitter_Density_NonOrtho (what a grid fit is reached through)
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
};

//! \brief The Dunlap fit on a Gaussian (non-orthogonal, Coulomb-metric) fit basis \a fbs.  THROWS when \a fbs has no Coulomb metric.
std::unique_ptr<DunlapDensityFit<double>> MakeDunlapDensityFitter(std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs);

//! \brief A density fit on a mesh: the Gaussian fit's face, plus the energy, which only its \c CoulombOperator can answer (V1.1b).
//! Its own class here, beside the hub.  Data-free.
class GridDensityFit
    : public virtual FunctionFitter_Density_NonOrtho<double>
    , public virtual FirstOrderDensityFit<double>
{
public:
    using FunctionFitter_Density_NonOrtho<double>::DoFit;   // the old face (3c retires it)
    using FirstOrderDensityFit<double>::DoFit;
    using FunctionFitter_Density_NonOrtho<double>::Repulsion;
    using FirstOrderDensityFit<double>::Repulsion;
    virtual double Integral() const override=0;            // one name for the two bases' Integral
};

//! \brief A density fit on \a deltaBasis (an orthogonal, \f$\delta\f$-like fit basis that answers its own Coulomb integrals,
//! \c BasisSet::PointCoulombIntegrals: a \c CoulombDeltaFit_IBS or a \c MultipoleFit_IBS).  THROWS when the basis is not orthogonal or
//! carries no Coulomb integrals.
std::unique_ptr<GridDensityFit>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis);

//! \brief A density fit on a PLANE-WAVE fit basis \a pw (orthonormal {G}): the fit is \f$\tilde\rho(G)\f$ itself, held privately;
//! the physics is the basis's Coulomb integrals.  THROWS when \a pw is not orthonormal or lacks the G-space faces.
std::unique_ptr<FirstOrderDensityFit<dcmplx>> MakeGSpaceDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS> pw);

} //namespace
