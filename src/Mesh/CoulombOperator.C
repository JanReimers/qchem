// File: Mesh/CoulombOperator.C  The Coulomb integrals of a point set: how the potential of charges on a mesh is found.
//
// MOVED HERE from qcFitting (V-ONE-OF-N 3b-ii, user ruling 2026-11-08 (f), doc/Records/VOneOfN_D2Design.md §1a): 1/r12 on a δ
// basis is that basis's COULOMB INTEGRALS, the physics a fitted Hartree term asks for -- not a fitter's private strategy.  The δ
// fit basis (qcBasisSet) owns one and answers for it (BasisSet::PointCoulombIntegrals), so the operator lives below qcBasisSet,
// beside the meshes and site multipoles it runs over.
//
// A fit basis says what a function is; a CoulombOperator says how the potential of a density on it is found (user ruling
// 2026-11-06, on lane S's design note).  The δ-fit routes (ρ fitted on a δ fit basis, whose points the XC mesh supplies) all end the
// same way -- v at the δ points, then V_ij = Σ_g c_g (ij|δ_g), the fit's coefficients contracted with the δ basis's exact projection
// integrals -- and differ only in how v is found, so that is the one thing that varies.
//
// MEASURED (lane S prototype, TeF₆ at a fixed density, doc/Records/ on the row): at a mesh's own rung the J error of any
// route that forms v at points is floored by the final contraction, 4.7e-4 per element at 80/29.  Analytic (ij|a) has no
// such floor.  The direct sum sits well above it (about 1e-2, flat in the radial count).
module;
#include <memory>
export module qchem.Mesh.CoulombOperator;
export import qchem.Mesh;                     // qcMesh::Mesh (the points a direct sum runs over)
export import qchem.Mesh.SiteMultipoles;      // what a radial solve runs over
import qchem.Mesh.Folded;                     // FoldedMesh: the orbits a folded operator scatters over
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
    //! \name THE DERIVATIVE of \c SelfRepulsion beside it, AT FIXED CHARGES (lane S 2026-11-23, the δ Hartree force's operator leg):
    //! \f$\partial(q^TKq)/\partial R_A\f$ per site of an atom-centred mesh, with \a q held and every point riding its owner site, so only the
    //! operator \f$K(R)\f$ moves -- scaled as the value is (no ½).  \f$q^TdKq\f$ is the exact derivative of the quadratic form for any K,
    //! symmetric or not: it sees only K's symmetric part.  The charges' own motion (the weights, the density at the points) is the
    //! caller's.  \c Can is false exactly when \c Make would throw.  Default: no pullback (the direct sum: its self term reads the weights;
    //! the multipole: its axis is channels, not points).  CLIENT: the grid density fitter's \c MakeHartreePullback.
    //!@{
    virtual bool       CanPullbackSelfRepulsion() const {return false;}
    virtual rvec3vec_t MakeSelfRepulsionPullback(const rvec_t& q) const;
    //!@}
};

//! \brief The direct sum over \a mesh's points: \f$v(r_g)=\sum_{h\ne g}q_h/|r_g-r_h|\f$ plus a self term.  \f$O(N^2)\f$,
//! matrix-free.  Copies the points: an operator outlives no borrowed mesh.
std::unique_ptr<CoulombOperator> MakeDirectSumCoulomb(const qcMesh::Mesh& mesh);

//! \brief The radial Poisson solve on an atom-centred mesh, to the \a sites' own lMax, for charges at the points.  The sites' radial
//! rule must be uniform in its index.  \f$O(N_{pts}N_{sites}(l+1)^2)\f$.
std::unique_ptr<CoulombOperator> MakeRadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites);

//! \brief The same solve for a multipole fit over \a sites: \a metric is that basis's diagonal.
std::unique_ptr<CoulombOperator> MakeMultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric);

//! \brief A point operator on the REPRESENTATIVES of a folded mesh: \a full runs over every point of \a folded's mesh, and this one answers
//! for one point per orbit (lane H 2026-11-24, the δ-fit Hartree on the representatives; OpenWork §3).
//!
//! WHY IT IS EXACT.  On an invariant mesh an operator built from it commutes with the group's action (the radial Poisson solve: its moments,
//! solve and synthesis, and their transposes, are covariant, complete in m), so an INVARIANT charge has an invariant potential.  The charge
//! at the representatives, \f$Q_o=W_o\rho_o\f$, then fixes the charge everywhere, \f$q_g=w_g\rho_o\f$ for each member g of orbit o, and the
//! potential is read back at the representative.  The energy \f$\tfrac12\sum_gq_gv_g=\tfrac12\sum_oQ_ov_o\f$ follows, as does the
//! self-repulsion pullback (taken on the full charges).  The CALLER owns the invariance of the charge: a SALC-imposed density on the mesh's
//! own group.  The cost is the full operator's: only the points the CLIENTS sample (collocation, contraction) fold.
//!
//! The orbit bookkeeping lives HERE, in the operator's wiring, so a fitter on the representatives sees an ordinary point basis whose
//! multiplicity is in its weights (the contract with lane S's analytic (ij|δ) fitter, 2026-11-24).  THROWS when \a folded has no fold, or
//! \a full is null.  \a full must be built on \a folded's mesh.
std::unique_ptr<CoulombOperator> MakeFoldedCoulomb(std::unique_ptr<CoulombOperator> full, const FoldedMesh& folded);

} //namespace
