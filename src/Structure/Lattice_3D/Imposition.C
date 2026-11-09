// File: Structure/Lattice_3D/Imposition.C  THE RUN'S SYMMETRY IMPOSITION: what a crystal run imposes, decided once from the
// lattice and the options (V-BS-SYMMETRY, user + controller 2026-11-02).
//
// DETECTION is the lattice's (Lattice_3D::GetSpaceGroup / ShubnikovOps); the DECISION -- which of the detected ops a run
// imposes, and on what (the k fold, the density star-average, the raster, the XC mesh, the spin pair) -- was resolved
// inside the GPW basis factory (DetectPointOps), and the plane-wave factory resolved its Shubnikov group a second time.  It
// lives here now, beside detection, so a basis is HANDED the ops it is built under instead of deciding them.
module;
#include <memory>
#include <vector>
export module qchem.Lattice_3D.Imposition;
export import qchem.Lattice_3D;
export import qchem.Symmetry.Lattice_3D.SpaceGroup;   // SymmetryPolicy, ReciprocalOp, DirectOp, SymOp
export import qchem.Symmetry.Imposition;              // SymmetryImposition: the face its consumers ask (S3)
export import qchem.Mesh.Folded;                      // qcMesh::FoldedMesh, MeshParams: what XCMeshSource makes
import qchem.Types;

export namespace qchem
{

//! What the caller asserted about the run's symmetry: the inputs \c ResolveImposition decides from.
struct ImpositionOptions
{
    bool             imposeSymmetry    = false;   //!< impose the detected group (impose-on-assert, doc/SymmetryUpgradePlan.md §3)
    std::vector<int> siteSpins         = {};      //!< collinear magnetic decoration per site (+1/-1/0); a staggering imposes Shubnikov
    bool             freeMeshInvariant = true;    //!< FREE run only: the XC mesh invariant under the detected group (KAxisPlan ruling 1)
    rvec3_t          kShift            = rvec3_t(0,0,0);   //!< the MP offset: the density group is restricted to ops mapping the mesh onto itself
    bool operator==(const ImpositionOptions&) const = default;
};

//! \brief The crystal point ops a run IMPOSES, from ONE space-group detection (or {} when not imposing -- the k fold and the
//! density star-average are one package, so a free run carries no density ops = trivial {E}).
//!
//! Faces and their consumers:
//!  - \c recipFold: reciprocal LINEAR ops U (+time reversal, τ=0 W-only guard) -- the basis folds its k-mesh to the IBZ.
//!    Time reversal supplies the inversion a non-symmorphic Td subgroup lacks, so the fold is maximal without the glide.
//!  - \c densityGroup: the full crystal point group with fractional translations {W|τ}, as the GROUP the composite density
//!    star-averages under (ρ̃ in G-space with the glide phase, the raster ρ with the voxel shift N·τ).  No time reversal.
//!  - \c directDensity {W|τ}: the same ops as a list, still read by the GPW blocks (the T1 {G} fold, the stream fold, the
//!    grid report) -- the remainder of V-SYM-SEQUESTER-OPS.
//!  - \c magneticDirect {W|τ,σ}: the Shubnikov group on a magnetic imposition ({} otherwise).  The two density faces are then
//!    its σ=None subgroup only (a sublattice swap is not a symmetry of ONE channel); σ is consumed by the spin-pair machinery.
//!  - \c meshInvariance: FREE run only, the detected spatial group the XC mesh is built invariant under (no density average).
//!  - \c spinPairOps: the Shubnikov group again, but only where every op maps every orbital block onto ITSELF -- a magnetic
//!    imposition of a single Γ block with no k shift -- so the (D↑, D↓) pair projection is defined; {} otherwise.
//! \brief An imposition that BUILDS the run's XC integration mesh under its own ops (V-BS-SYMMETRY S2, controller 2026-11-03).
//! A BEHAVIOUR face, not a data one: the XC quadrature builder asks for the finished mesh, and no op set or cell matrix leaves
//! the imposition.  The builder reaches it by a cross-cast from \c Symmetry::SymmetryImposition; an imposition that does not
//! realize it (none yet, a molecule's point group later) leaves the builder on the Structure's own plain mesh.
class XCMeshSource
{
public:
    virtual ~XCMeshSource() = default;
    //! The XC integration mesh at \a mp: invariant under the imposed ops and orbit-folded (plus the Shubnikov σ tags and the
    //! flip-fixed flags on a magnetic imposition) when the run imposes; invariant but UNFOLDED on a free run with mesh
    //! invariance (the density is never averaged); the cell's plain mesh otherwise.
    virtual qcMesh::FoldedMesh MakeXCMesh(const qcMesh::MeshParams& mp) const=0;
};

//! For a SYMMORPHIC crystal every τ=0 and the density faces collapse to the W-only ops.
//! It REALIZES \c Symmetry::SymmetryImposition, the face the composite density and wave function ask (V-BS-SYMMETRY S3), and
//! \c XCMeshSource, the face the XC quadrature builder asks (S2).
class CrystalImposition : public virtual Symmetry::SymmetryImposition, public virtual XCMeshSource
{
public:
    const std::vector<Symmetry::SymOp>&               SpinPairOps   () const override {return spinPairOps;}
    std::shared_ptr<const Symmetry::StarGroup> DensityGroup() const override {return densityGroup;}
    qcMesh::FoldedMesh MakeXCMesh(const qcMesh::MeshParams& mp) const override;

    ImpositionOptions                                    options;        //!< what this was resolved from (a consumer handed both checks they agree)
    Symmetry::Lattice_3D::SymmetryPolicy                 policy;         //!< the ONE run-level imposition decision (§3)
    std::vector<Matrix3D<double>>                        recipFold;      //!< linear reciprocal ops (+TR) for the k-fold ({} unless imposed, grey)
    std::shared_ptr<const Symmetry::StarGroup>           densityGroup;   //!< the density's star-average group (null = {E}: not imposed)
    std::vector<Symmetry::Lattice_3D::DirectOp>          directDensity;  //!< the same {W|τ} as a list, for the GPW blocks ({} unless imposed)
    std::vector<Symmetry::Lattice_3D::SymOp>             magneticDirect; //!< the Shubnikov group of a magnetic imposition ({} = grey or free)
    std::vector<Symmetry::Lattice_3D::SymOp>             meshInvariance; //!< free run: the XC mesh's invariance group ({} when imposing)
    std::vector<Symmetry::SymOp>                         spinPairOps;    //!< magneticDirect on a single unshifted Γ block, else {}
    //! The cell the ops were resolved on, which \c MakeXCMesh builds the mesh over.  A COPY, held so the imposition owns its
    //! own lifetime: it outlives the \c Lattice_3D it was resolved from inside every basis that holds it.
    std::shared_ptr<const UnitCell>                      cell;
};

//! \brief THE ONE PLACE a crystal run decides what it imposes (doc/SymmetryUpgradePlan.md §3).  Detection always runs (it is
//! the lattice's, and cached there); the imposed faces fill only when \a o.imposeSymmetry.  A genuine staggering in
//! \a o.siteSpins resolves to the magnetic (Shubnikov) group -- the grey group would star-average +m onto -m sites and erase
//! the order.  Every imposed set passes the IMPOSITION GATE (closure: a star-average is a projector only over a group), so an
//! open set throws here rather than diverging later (the Si 2x2x2 supercell, 2026-09-07).
CrystalImposition ResolveImposition(const Lattice_3D& lat, const ImpositionOptions& o);

} // namespace qchem
