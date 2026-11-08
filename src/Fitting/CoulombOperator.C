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

export namespace qchem::Fitting
{

//! \brief The Coulomb potential of a density on a fit basis.  Data-free.  \a q and the answer are both indexed like the fit basis's
//! functions.  \a q is \f$\langle f_a|\rho\rangle\f$, so on a \f$\delta\f$ basis \f$q_g=w_g\rho(r_g)\f$, a CHARGE, and an operator
//! never needs the density's own values.  The answer is the potential's coefficients on the same functions: its values at the points.
class CoulombOperator
{
public:
    virtual ~CoulombOperator() = default;
    virtual rvec_t Potential(const rvec_t& q) const=0;   //!< the coefficients of \f$v\f$
    //! \f$\langle\tilde\rho|\tilde\rho\rangle\f$ for \a q and its own \a v.  On points that is \f$\sum_gq_gv_g\f$; an operator whose
    //! functions overlap between sites owes its own answer.
    virtual double SelfRepulsion(const rvec_t& q, const rvec_t& v) const
    {
        double e=0.0;
        for (size_t g=0; g<q.size(); g++) e+=q[g]*v[g];
        return e;
    }
    //! \brief The Hartree energy of \a q and its own \a v.  \a trDV is \f$\mathrm{Tr}(DV)\f$, which only the density can supply.  On points
    //! both legs of \f$\mathrm{Tr}(DV)-\tfrac12\langle\tilde\rho|\tilde\rho\rangle\f$ are the same number, so it is \f$\tfrac12\sum_gq_gv_g\f$
    //! and \a trDV is not needed: first order in the mesh (ruled 2026-11-06).
    virtual double HartreeEnergy(double /*trDV*/, const rvec_t& q, const rvec_t& v) const {return 0.5*SelfRepulsion(q,v);}
    //! Does \c HartreeEnergy read \a trDV?  Only an operator whose two legs differ does (the multipole); the caller then
    //! pays for \f$\mathrm{Tr}(DV)\f$, a full matrix contraction, only when it is used.
    virtual bool UsesTrDV() const {return false;}
};

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

//! \brief The direct sum over \a mesh's points: \f$v(r_g)=\sum_{h\ne g}q_h/|r_g-r_h|\f$ plus a self term.  \f$O(N^2)\f$,
//! matrix-free.  Copies the points: an operator outlives no borrowed mesh.
std::unique_ptr<CoulombOperator> MakeDirectSumCoulomb(const qcMesh::Mesh& mesh);

//! \brief The radial Poisson solve on an atom-centred mesh, to the \a sites' own lMax, for charges at the points.  The sites' radial
//! rule must be uniform in its index.  \f$O(N_{pts}N_{sites}(l+1)^2)\f$.
std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites);

//! \brief The same solve for a multipole fit over \a sites: \a metric is that basis's diagonal.
std::unique_ptr<CoulombOperator> MakeMultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric);

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

//! \brief A density fit on \a deltaBasis (an orthogonal, \f$\delta\f$-like fit basis), with \a op finding the potential.
//! THROWS when \a deltaBasis is not orthogonal: a diagonal metric is what makes the projection the fit.
std::unique_ptr<GridDensityFit>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis, std::unique_ptr<CoulombOperator> op);

//! \brief A density fit on a PLANE-WAVE fit basis \a pw (orthonormal {G}): the fit is \f$\tilde\rho(G)\f$ itself, held privately;
//! the physics is the basis's Coulomb integrals.  THROWS when \a pw is not orthonormal or lacks the G-space faces.
std::unique_ptr<FirstOrderDensityFit<dcmplx>> MakeGSpaceDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS> pw);

} //namespace
