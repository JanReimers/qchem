// File: BasisSet/Fit_Types.C  The VOCABULARY of the fit-basis factories -- not a fit-basis interface.
//
// Two plain value types that a factory takes in or hands back.  They lived in Fit_IBS.C until 2026-08-23,
// which was wrong for a simple reason (user): that file DEFINES THE FIT-BASIS INTERFACES, and neither of
// these is one -- no fit-basis face names either of them.  What they belong to is the CREATION story:
//
//   VxcFit        -- goes IN:  which representation should be built?
//   FitQuadrature -- comes OUT: the finished real-space quadrature a delta basis is built over.
//
// A LEAF module, deliberately: it needs the mesh and the symmetry fold and nothing else, so it can sit
// below Orbital_DFT_IBS and tBasisSet (whose fit factories take it in) without either importing the fit-basis faces to get
// at it.
module;
#include <memory>
#include <vector>   // FitQuadrature::sigmas/flipFixed (Shubnikov S3)
#include <string>
export module qchem.BasisSet.Fit_Types;
export import qchem.Mesh;                       // qcMesh::Mesh / MeshParams
export import qchem.Mesh.Folded;                // qcMesh::FoldedMesh -- WHAT FitQuadrature now IS
export import qchem.Symmetry.Lattice_3D.Fold;   // Fold -- re-exported: callers still build one to hand over
import qchem.Symmetry.Group;                    // Frame / GroupPtr: what a fit family FORCES (Forces, Group face step 3b)

export namespace qchem::BasisSet
{

//! \brief A FINISHED real-space XC quadrature (doc/SymmetryUpgradePlan.md §6a): the mesh, plus its orbit
//! \c fold when the run imposes symmetry (empty = free run).  This is what the DELTA-fit XC route consumes
//! -- a pure quadrature, deliberately NOT a fit basis.
//! Produced by the run point's \c qcMesh::MeshMaker (V-MESH-FACTORY §4.1): the low-level mesh work (grid build, invariance under the
//! imposed ops, fold preparation) is the maker's, over the geometry it owns, with the ops from the run's \c qcMesh::SymmetrySource -- not the basis's, and not
//! the Hamiltonian assembly's, which only hands the result to the δ basis and the XC term.
//!
//! ★ **IT IS NOW `qcMesh::FoldedMesh` (2026-09-08), and the alias is the point.**  This was four loose
//! members -- mesh, fold, sigmas, flipFixed -- declared in a FITTING header, wired together by hand at
//! every use, with nothing checking that the fold belonged to the mesh.  But the pointwise star-average is
//! an EXACT projector only because the mesh is invariant under the ops the fold came from: that is a joint
//! property of the two, so they are one type, the pairing is checked once in its constructor, and the
//! star-average is a member instead of something each consumer re-writes (user, 2026-09-08).
//! The name stays as the FIT-FACTORY's vocabulary for "the quadrature you asked me to build"; the TYPE is
//! the mesh library's.
using FitQuadrature = qcMesh::FoldedMesh;

//! \brief WHICH REPRESENTATION carries \f$v_{xc}\f$ -- the argument that picks between the fit-basis types,
//! and therefore an argument of \c CreateVxcFitBasisSet.
//!  - \c Delta: the \f$\delta\f$ basis (\c DeltaFit_IBS) -- \c n_pts delta functions on the quadrature
//!    mesh, diagonal metric, fit = identity, so the "coefficients" ARE the point values.
//!  - \c PlaneWave: plane waves on the \f$\{G\}\f$ ball -- the PERIODIC lineage's fitted representation (band-limited; the
//!    projection is an FFT on its raster).
//!  - \c Auxiliary: the non-periodic lineage's own fitted AUXILIARY basis -- the Gaussian A1 set on a molecule, the radial
//!    fit set on an atom (2026-10-31: named when δ became the molecular / atom default, so the fit stays selectable -- user
//!    ruling "δ default, A1 kept", and pin 30: it keeps the framework honest over EVERY fit basis).
//!  - \c Auto: the caller did not choose.  Resolved by the POLICY layer -- the Hamiltonian: the periodic one picks Delta or
//!    PlaneWave (it alone knows whether the run is polarized / on a Becke grid), the molecular ones pick Delta (the user's
//!    default, 2026-10-31).  If it reaches a factory unresolved it is read as "your own fitted representation" (the
//!    historical pairing; \c qcMesh::UnitCellKind::Auto falls through the same way).
//!  - \c Multipole: radial channels per (site, l, m) on the same mesh (\c MultipoleFit_IBS).  Finite structures only, and the
//!    Hamiltonian builds it over a \f$\delta\f$ basis; no basis factory reads this value.
//! A combination a lineage cannot realize (a plane-wave fit on a molecule, an auxiliary fit on a cell) THROWS at the factory,
//! naming it: an open gap of the product space (pin 30), never a silent substitute.
//!
//! It does not live in \c qcMesh (where it briefly did, 2026-08-22): \c MeshParams describes POINTS, and
//! this describes FUNCTIONS -- the two are exactly the orthogonal axes the XC separation is about, so
//! folding one into the other's parameter block would have re-welded them at the type level while the code
//! was busy unwelding them.  It does not live in qcHamiltonian either, because the factory that reads it is
//! here.  \c qchem::Hamiltonian::VxcFit is an alias, so the user-facing spelling is unchanged.
enum class VxcFit {Auto, PlaneWave, Delta, Auxiliary, Multipole};

//! \brief WHAT A Vxc FIT FAMILY FORCES on the potential it carries, as a Group in \a frame (doc/Pins.md pin 30; Group face step 3b).  The
//! family fixes the group's TYPE and the frame supplies where -- nothing is built:
//!  - the plane-wave ball forces T of the frame's lattice (and cannot exist without one);
//!  - the lineage's A1 set forces the frame's own family: lattice-summed on a cell (T of its lattice), the l = 0 set about the atom's site
//!    (O(3) there), atom-centred Gaussians of every l on a molecule ({e}).  The lineage follows the structure's kind, so the frame its group
//!    offers (\c Symmetry::Group::MakeFrame) names it;
//!  - \f$\delta\f$ functions, site multipoles and an unresolved \c Auto force nothing ({e}).
//! It FAILS, with the reason, when the frame lacks what the family needs.  ONE source of truth: every fit class answers its own
//! \c FIT_SF_ABS::Forces from this, in its own frame, and the pin-30 pre-flight (\c Hamiltonian::AdmitVxc) asks it of the structure's
//! group's frame.
Outcome<Symmetry::GroupPtr,std::string> Forces(VxcFit fit, const Symmetry::Frame& frame);

}//namespace
