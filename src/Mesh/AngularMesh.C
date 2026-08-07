// File: AngularMesh.C  Angular quadrature mesh (concrete value type) + typed factory.
//
// Unit directions Omega_i and weights w_i normalised so that  sum_i w_i = 4*pi.  The schemes
// (Lebedev / GaussLegendre) are plain builder FUNCTIONS that return an AngularMesh --
// no class hierarchy.
module;
#include <utility>
#include <vector>
export module qchem.Mesh.Angular;
export import qchem.Types;
export import qchem.Mesh;            // AngularKind, MeshParams

export namespace qchem::qcMesh
{

class AngularMesh
{
public:
    AngularMesh() = default;
    AngularMesh(rvec3vec_t d, rvec_t w) : itsD(std::move(d)), itsW(std::move(w)) {}
    const rvec3vec_t& Dirs() const {return itsD;}   //!< unit directions Omega_i
    const rvec_t&     W   () const {return itsW;}   //!< weights, sum = 4*pi
    size_t          size  () const {return itsD.size();}
private:
    rvec3vec_t itsD;
    rvec_t     itsW;
};

//! \brief Build an angular mesh of the requested kind from the typed parameters.
AngularMesh MakeAngular(const MeshParams&);

//! \brief One tabulated Lebedev rule: its direction count and the polynomial degree it GUARANTEES.
//!
//! \a degree is the CONSTRUCTED degree, not the measured one.  Two rules over-deliver when scanned
//! monomial-by-monomial (302 measures 31, 434 measures >=40), because monomials odd under the rule's
//! octahedral symmetry vanish identically on both sides and never discriminate.  Guaranteeing less than
//! you deliver is honest; the reverse is not -- so the table states the guarantee.
struct LebedevRule {int nDir; int degree;};

//! \brief The Lebedev menu, ascending in degree.  Public because the degree ladder has PRINCIPLED GAPS
//! that a caller may need to see: degrees 9, 13, 25 and 27 are absent because those orders carry NEGATIVE
//! WEIGHTS (74, 230, 266 -- excluded by the generator audit) or shipped a weight-sum bug (the 32-point
//! degree-9 rule, removed).  So the jump from degree 23 to 29 costs 194 -> 302 directions, +55%.
const std::vector<LebedevRule>& LebedevMenu();

//! \brief The cheapest Lebedev rule delivering AT LEAST \a degree.  Rounds UP -- the only rule that
//! guarantees the requested exactness -- and ANNOUNCES the substitution when it has to round, naming the
//! gap, so a caller never silently receives a 55%-dearer grid than it asked for.  Throws past the menu.
LebedevRule ResolveLebedev(int degree);

//! \brief Which HIGH-SYMMETRY directions a rule places points on: the \f$\langle100\rangle\f$ axes,
//! the \f$\langle110\rangle\f$ edges, the \f$\langle111\rangle\f$ body diagonals.
//!
//! Two rules can share a degree and still differ HERE, and that difference is load-bearing: a quadrature
//! direction lying along a bond is exactly what \c MeshParams::angRot exists to steer away from (§6a).
//! The 6- and 8-direction rules are the lowest-order instance (both degree 3, on \f$\langle100\rangle\f$
//! and \f$\langle111\rangle\f$ respectively) and higher-order ties behave the same way -- so this is
//! MEASURED from the directions rather than annotated per rule.  This file's history is the argument for
//! computing it: the hand-written degree annotations had drifted on two rules before they were measured.
//!
//! \warning The frame is CUBIC.  These are the octahedral classes, so "none of the three" means GENERIC
//! WITH RESPECT TO A CUBIC LATTICE -- not structureless.  The 12-direction rule is the case in point: it
//! is the ICOSAHEDRON (\f$r=\varphi/\sqrt{1+\varphi^2}\f$, \f$s=1/\sqrt{1+\varphi^2}\f$), which has
//! its own special directions (12 vertices, 20 face centres, 30 edge midpoints) that this classifier
//! cannot see.  For the site-adapted work that is the RIGHT reading -- a rule generic w.r.t. the lattice
//! cannot accidentally put a quadrature point on a bond -- but do not read "none" as "no symmetry".
struct SpecialOrbits {bool axes100=false, edges110=false, corners111=false;};
SpecialOrbits ClassifyOrbits(const AngularMesh&, double tol=1e-9);

//! \brief Build ONE tabulated Lebedev rule, addressed by its direction count.
//!
//! The interface speaks DEGREES (\c MeshParams::angularDegree) and the builder speaks RULES;
//! \c ResolveLebedev is the bridge.  Exported because a rule is a real object worth addressing
//! directly -- the transcription audits check each tabulated rule, not the degree mapping.
AngularMesh LebedevAngular(int numDir);

} //export namespace qchem::qcMesh

// Per-scheme builders -- declared NON-exported here, implemented in separate files; only MakeAngular
// (above) calls them.
namespace qchem::qcMesh
{
    AngularMesh GaussLegendreAngular(int L);
}
