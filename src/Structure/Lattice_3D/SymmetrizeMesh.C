// File: Structure/Lattice_3D/SymmetrizeMesh.C  Symmetry folding of quadrature meshes (T2 groundwork).
//
// doc/SymmetryUpgradePlan.md §2a: qcSymmetry owns the FOLD (a pure partition over primitive vectors);
// the Mesh-typed convenience lives CLIENT-side.  It sits HERE, in qcStructure, where Mesh and Symmetry
// MEET WITHOUT DEPENDING ON EACH OTHER (user decision 2026-08-01): qcMesh stays a pure quadrature leaf,
// qcSymmetry stays the pure group-theory leaf, and qcStructure -- which links both, builds the periodic
// Becke mesh (UnitCell), and owns the crystal's detected SpaceGroup (Lattice_3D) -- is where quadrature
// meets group action.  (The G-space sibling SymmetrizeGMap lives with the G-currency in qcBasisSet.)
//
// The {r}-quadrature toolkit for the T2 site-adapted Becke grid, BEFORE any XC wiring:
//   - FoldMesh       : the partition itself (consumers needing per-member data, e.g. SymmetrizeValues);
//   - Symmetrize     : the compact folded mesh -- one representative per star, weight = the SUM of the
//                      member weights (exact for a totally symmetric integrand whatever the weight
//                      distribution, == w*n_star when the mesh is op-invariant);
//   - UnmatchedCounts: the INVARIANCE checker -- per op, how many points fail to land on a mesh point.
//                      A pointwise star-average / fold is exact only on an invariant mesh (§5/§6 T2);
//   - MakeInvariant  : the CONSTRUCTIVE half -- group-average any mesh into an invariant one (the
//                      union of all op images at weight w/|ops|, coincident points merged).  Each
//                      rotated copy is a quadrature of the same polynomial exactness, so the average
//                      keeps the original exactness degree AND gains exact op-invariance; folding it
//                      back (Symmetrize) returns to ~the original point count.  The generic route to
//                      the site-group-adapted Becke grid before any bespoke minimal-direction design.
//
// Mesh points are CARTESIAN (Mesh's currency); the ops act on FRACTIONAL coordinates on the unit
// torus, so every entry point takes the cell matrix A (columns = lattice vectors, r = A f) -- or a
// SpaceGroup, which carries both its ops and its cell (§2c: no ops at the call site).
//
// The interface holds declarations only; every body (and the module-private helpers they call) is in Imp/SymmetrizeMesh.C (CLAUDE.md:
// bodies in Imp/).  They were inline here, calling `static` helpers: an exported inline exposing an internal-linkage entity is
// ill-formed ([basic.link] exposure), and clang 23 is the first to say so (lane G, 2026-11-23).
module;
#include <vector>
export module qchem.SymmetrizeMesh;
export import qchem.Mesh;                            // qcMesh::Mesh
export import qchem.Symmetry.Lattice_3D.Fold;        // Fold, SymOp, FoldPointsPeriodic, SymmetrizeValues
export import qchem.Symmetry.Lattice_3D.SpaceGroup;  // SpaceGroup (the ops-encapsulated overloads)

export namespace qchem
{

namespace SL = qchem::Symmetry::Lattice_3D;

//! \brief Default coincidence tolerance for mesh-point matching, in FRACTIONAL coordinates (the torus
//! metric these routines work in), so it is cell-size independent by construction.
//!
//! Named once rather than repeated as a literal at each entry point: every routine here decides the SAME
//! question -- "are these two mesh points the same point?" -- so they must agree, and a reader should not
//! have to diff four default arguments to confirm that they do (R2.12).  1e-8 is far below any real mesh
//! spacing and far above the round-off of a fractional coordinate that has been through one W matrix.
inline constexpr double kMeshMatchTol = 1e-8;

//! The fold partition of \a m's points under \f$\{W|\tau\}\f$ ops acting on fractional coordinates
//! (torus metric, tolerance \a tol in fractional units).
SL::Fold FoldMesh(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                         const std::vector<SL::SymOp>& ops, double tol);

//! The compact folded mesh: representatives only, weight = the sum over the star's members --
//! exact for a totally symmetric integrand:
//! \f$\sum_i w_i f(r_i) = \sum_\mathrm{stars} f(r_\mathrm{rep})\sum_\mathrm{members} w\f$.
qcMesh::Mesh Symmetrize(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                               const std::vector<SL::SymOp>& ops, double tol);

//! Per-op invariance defect: how many of \a m's points an op maps OFF the mesh (all zeros ==
//! the mesh is invariant, the §6 T2 precondition for pointwise star-averaging and folding).
std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                                        const std::vector<SL::SymOp>& ops, double tol);

//! Group-average \a m into an op-invariant mesh: the union of all op images, each carrying weight
//! \f$w/|\mathrm{ops}|\f$, coincident points (torus metric, \a tol) merged.  Pass the FULL op set
//! including the identity -- the average is a projector only over a group.  Preserves \f$\sum w\f$
//! and the quadrature's polynomial exactness; the point count grows toward \f$|\mathrm{ops}|\times\f$
//! for a generic mesh (\c Symmetrize folds it back to ~the original count).
qcMesh::Mesh MakeInvariant(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                                  const std::vector<SL::SymOp>& ops, double tol);

// ---- The ops-encapsulated overloads (§2c): the SpaceGroup supplies its own {W|τ} set + cell. ----
SL::Fold         FoldMesh       (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);
qcMesh::Mesh     Symmetrize     (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);
std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);
qcMesh::Mesh     MakeInvariant  (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);

//! \brief SITE-GROUP-INVARIANT angular quadrature of polynomial degree \a L (doc/SymmetryUpgradePlan.md
//! §6a W2 -- "the site group fixes both how many and which directions").  Directions = unions of
//! orbits of directions under the site ops, so the grid is op-invariant BY CONSTRUCTION with equal
//! weights within each orbit.  THE MIXED RULE (§6a W2 refinement): the candidate pool holds the
//! group's SPECIAL orbits first (op axes + mirror-plane directions -- the small Lebedev-like orbits),
//! filtered against \a avoid (the ATOM's actual bond directions -- only some special families are
//! poison, e.g. diamond T_d's bonded <111> tetrahedron, while <100>/<110>/the anti-bond tetrahedron
//! point between bonds and are safe), then generic Fibonacci-sphere seeds.  Orbit weights solve the
//! spherical-harmonic moment conditions \f$Q(\bar Y_{lm}) = \sqrt{4\pi}\,\delta_{l0}\,\delta_{m0},\
//! l\le L\f$ by NNLS (Lawson-Hanson, per-DIRECTION cost-biased column selection), which guarantees
//! nonnegative weights and picks a SPARSE orbit support -- plain least squares was measured to leave
//! a ~-5e-4 negative-weight floor at L=29 however many orbits were added.  Only the positive-weight
//! support is emitted; \f$\sum w = 4\pi\f$ (the \f$Y_{00}\f$ condition).  Deterministic.
//! \param ops   the site group as CARTESIAN orthogonal matrices (rotations + reflections; include E).
//!              Fractional \f$W\f$ ops convert via \f$A W A^{-1}\f$ upstream.
//! \param L     the polynomial exactness degree (GL-\a L quality class).
//! \param avoid unit directions (e.g. the atom's bonds) no quadrature point may lie on.
qcMesh::Mesh MakeInvariantAngularMesh(const std::vector<Matrix3D<double>>& ops, int L,
                                      const std::vector<rvec3_t>& avoid = {});

} // export namespace
