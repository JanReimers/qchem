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
};

//! \brief The direct sum over \a mesh's points: \f$v(r_g)=\sum_{h\ne g}q_h/|r_g-r_h|\f$ plus a self term.  \f$O(N^2)\f$,
//! matrix-free.  Copies the points: an operator outlives no borrowed mesh.
std::unique_ptr<CoulombOperator> MakeDirectSumCoulomb(const qcMesh::Mesh& mesh);

//! \brief The radial Poisson solve on an atom-centred mesh, to the \a sites' own lMax, for charges at the points.  The sites' radial
//! rule must be uniform in its index.  \f$O(N_{pts}N_{sites}(l+1)^2)\f$.
std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites);

//! \brief The same solve for a multipole fit over \a sites: \a metric is that basis's diagonal.
std::unique_ptr<CoulombOperator> MakeMultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric);

//! \brief A density fit on \a deltaBasis (an orthogonal, \f$\delta\f$-like fit basis), with \a op finding the potential.
//! It answers the same face the Gaussian fit answers, so \c FittedCD and \c FittedVee cannot tell them apart.
//! THROWS when \a deltaBasis is not orthogonal: a diagonal metric is what makes the projection the fit.
std::unique_ptr<FunctionFitter_Density_NonOrtho<double>>
MakeGridDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis, std::unique_ptr<CoulombOperator> op);

} //namespace
