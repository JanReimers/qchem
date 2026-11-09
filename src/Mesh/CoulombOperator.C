// File: Mesh/CoulombOperator.C  The Coulomb integrals of a point set: how the potential of charges on a mesh is found.
//
// MOVED HERE from qcFitting (V-ONE-OF-N 3b-ii, user ruling 2026-11-08 (f), doc/Records/VOneOfN_D2Design.md §1a): 1/r12 on a δ
// basis is that basis's COULOMB INTEGRALS, the physics a fitted Hartree term asks for -- not a fitter's private strategy.  The δ
// fit basis (qcBasisSet) owns one and answers for it (BasisSet::PointCoulombIntegrals), so the operator lives below qcBasisSet,
// beside the meshes and site multipoles it runs over.
//
// A fit basis says what a function is; a CoulombOperator says how the potential of a density on it is found (user ruling
// 2026-11-06, on lane S's design note).  The grid routes (a density sampled on a mesh) all end the same way -- v at the mesh's
// points, then the same quadrature XC uses -- and differ only in how v is found, so that is the one thing that varies.
//
// MEASURED (lane S prototype, TeF₆ at a fixed density, doc/Records/ on the row): at a mesh's own rung the J error of any
// route that forms v at points is floored by the final contraction, 4.7e-4 per element at 80/29.  Analytic (ij|a) has no
// such floor.  The direct sum sits well above it (about 1e-2, flat in the radial count).
module;
#include <memory>
export module qchem.Mesh.CoulombOperator;
export import qchem.Mesh;                     // qcMesh::Mesh (the points a direct sum runs over)
export import qchem.Mesh.SiteMultipoles;      // what a radial solve runs over
export import qchem.Types;

export namespace qchem::qcMesh
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

//! \brief The direct sum over \a mesh's points: \f$v(r_g)=\sum_{h\ne g}q_h/|r_g-r_h|\f$ plus a self term.  \f$O(N^2)\f$,
//! matrix-free.  Copies the points: an operator outlives no borrowed mesh.
std::unique_ptr<CoulombOperator> MakeDirectSumCoulomb(const qcMesh::Mesh& mesh);

//! \brief The radial Poisson solve on an atom-centred mesh, to the \a sites' own lMax, for charges at the points.  The sites' radial
//! rule must be uniform in its index.  \f$O(N_{pts}N_{sites}(l+1)^2)\f$.
std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites);

//! \brief The same solve for a multipole fit over \a sites: \a metric is that basis's diagonal.
std::unique_ptr<CoulombOperator> MakeMultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric);

} //namespace
