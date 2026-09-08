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
// below Orbital_DFT_IBS and tBasisSet (both of which declare CreateXCQuadrature) without either importing
// the fit-basis faces to get at it.
module;
#include <memory>
#include <vector>   // FitQuadrature::sigmas/flipFixed (Shubnikov S3)
export module qchem.BasisSet.Fit_Types;
export import qchem.Mesh;                       // qcMesh::Mesh / MeshParams
export import qchem.Mesh.Folded;                // qcMesh::FoldedMesh -- WHAT FitQuadrature now IS
export import qchem.Symmetry.Lattice_3D.Fold;   // Fold -- re-exported: callers still build one to hand over

export namespace qchem::BasisSet
{

//! \brief A FINISHED real-space XC quadrature (doc/SymmetryUpgradePlan.md §6a): the mesh, plus its orbit
//! \c fold when the run imposes symmetry (empty = free run).  This is what the DELTA-fit XC route consumes
//! -- a pure quadrature, deliberately NOT a fit basis.
//! Produced by \c CreateXCQuadrature so all the low-level mesh work (grid build, group-averaging it
//! invariant under the imposed ops, fold preparation) lives with the basis factories -- which own the
//! cell and the §3 policy -- not in the Hamiltonian assembly.
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
//!  - \c PlaneWave: the lineage's own FITTED representation -- plane waves on the \f$\{G\}\f$ ball for a
//!    periodic basis (band-limited; the projection is an FFT on its raster).
//!  - \c Auto: the caller did not choose.  Resolved by the POLICY layer -- the Hamiltonian, which is the
//!    only layer that knows whether the run is polarized -- and, if it reaches a factory unresolved, read
//!    as "your own fitted representation" (the historical pairing; \c qcMesh::UnitCellKind::Auto falls
//!    through the same way).  A molecular basis, which has only its Gaussian auxiliary representation,
//!    is asked with \c Auto and never needs to interpret the other two.
//!
//! It does not live in \c qcMesh (where it briefly did, 2026-08-22): \c MeshParams describes POINTS, and
//! this describes FUNCTIONS -- the two are exactly the orthogonal axes the XC separation is about, so
//! folding one into the other's parameter block would have re-welded them at the type level while the code
//! was busy unwelding them.  It does not live in qcHamiltonian either, because the factory that reads it is
//! here.  \c qchem::Hamiltonian::VxcFit is an alias, so the user-facing spelling is unchanged.
enum class VxcFit {Auto, PlaneWave, Delta};

}//namespace
