// File: Mesh/Maker.C  qcMesh::MeshMaker -- the abstract seam "build me a mesh from this recipe" (V-MESH-FACTORY, doc/MeshFactoryPlan.md §3).
//
// A library that must CHOOSE a recipe itself (the GPW basis's PP-fallback mesh at its own projector reach, the Hamiltonian factory's
// PP quadrature, the SAD seed) holds one of these, and nothing else of the structure's mesh machinery: it names only a recipe and a
// built mesh, so the holder imports neither the structure factory nor a concrete structure.  The concrete realization is
// qchem.Structure.MeshFactory's StructureMeshMaker<S>, built at the edge (the setup) where the structure's type is known -- the
// house dependency-inversion rule (CLAUDE.md, "Constructors pass concrete objects; factories concentrate the concrete associations").
//
// WHY A MODULE OF ITS OWN and not qchem.Mesh.Recipe (where §3 first put it): Make returns a qcMesh::Mesh, and qchem.Mesh imports the
// recipe, so the seam has to sit above both.
module;
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
export module qchem.Mesh.Maker;
export import qchem.Mesh;          // qcMesh::Mesh, the built object Make returns
export import qchem.Mesh.Recipe;   // qcMesh::MeshParams, the recipe Make reads
export import qchem.Mesh.Folded;   // qcMesh::FoldedMesh, what the symmetric Make returns (BasisSet::FitQuadrature)
export import qchem.Symmetry.SymOp;   // Symmetry::SymOp, what a SymmetrySource hands the maker

export namespace qchem::qcMesh
{

//! \brief WHAT A MESH MUST RESPECT: the ops it is built invariant under, and whether the run averages the density over their orbits
//! (controller + user rulings 2026-11-22, doc/MeshFactoryPlan.md §4.1 as built).  DATA-FREE and geometry-free: a source supplies OPS only,
//! the \c MeshMaker owns the geometry they act on, so a source resolved at one cell can never move a mesh built over another.  Realized by
//! \c CrystalImposition (a crystal run), \c NoSymmetry (nothing asked), and \c MolecularPointGroupSource (a symmetric molecule: its D2h-type subgroup, OOD-SOLID V-FOLD-NEUTRAL).
class SymmetrySource
{
public:
    virtual ~SymmetrySource()=default;
    //! The ops the mesh must be invariant under ({} = none: the structure's plain mesh).
    virtual std::vector<Symmetry::SymOp> MeshOps() const=0;
    //! Does the run star-average the density over the ops' orbits?  Yes: the maker FOLDS the invariant mesh (the orbit partition the
    //! averaging consumes).  No: an invariant mesh for symmetric numerics, unfolded, and the density is never averaged (a free run).
    virtual bool AveragesDensityOverOrbits() const=0;
};

//! The null source: no ops, no averaging.  A molecule's, an atom's, the plane-wave family's.
class NoSymmetry : public SymmetrySource
{
public:
    std::vector<Symmetry::SymOp> MeshOps() const override {return {};}
    bool AveragesDensityOverOrbits() const override {return false;}
};

//! \brief Builds a real-space integration mesh from a recipe over ONE fixed structure (the one the realization closes over).
//! Two overloads (the user, 2026-11-22: one interface, not two parallel class networks for raw and folded meshes): the plain mesh, and
//! the mesh a \c SymmetrySource asks for.  The second is a TEMPLATE METHOD: the branch logic is the same for every structure; only the
//! two steps that need the geometry -- the invariant mesh, and the invariant mesh with its orbit fold (a cell's matrix, its Shubnikov
//! flip-fixed points) -- are the realization's protected hooks.  Neither hook is reached when the source has no ops.
class MeshMaker
{
public:
    virtual ~MeshMaker()=default;
    //! The mesh of the bound structure at recipe \a mp.  Throws what the structure's builder throws (a pin-30 refusal, a broken invariant).
    virtual Mesh Make(const MeshParams& mp) const=0;
    //! The mesh at \a mp that \a sym asks for: plain when it has no ops; invariant under them, unfolded, when the run does not average
    //! the density; invariant AND orbit-folded when it does.
    FoldedMesh Make(const MeshParams& mp, const SymmetrySource& sym) const
    {
        const std::vector<Symmetry::SymOp> ops=sym.MeshOps();
        if (ops.empty()) return {std::make_shared<const Mesh>(Make(mp)), {}};
        if (sym.AveragesDensityOverOrbits()) return MakeFolded(mp, ops);
        Mesh mesh=MakeInvariant(mp, ops);
        // FREE run: the same invariant mesh an imposed run gets, but NO fold -- nothing averages the density.
        std::cout << "[XC mesh] free run: mesh invariant under the " << ops.size() << " detected ops, "
                  << mesh.size() << " points; the density is NOT averaged" << std::endl;
        return {std::make_shared<const Mesh>(std::move(mesh)), {}};
    }
protected:
    //! The bound structure's mesh at \a mp, invariant under \a ops (non-empty).
    virtual Mesh       MakeInvariant(const MeshParams& mp, const std::vector<Symmetry::SymOp>& ops) const=0;
    //! The bound structure's mesh at \a mp, invariant under \a ops (non-empty), WITH its orbit partition (and on a magnetic op set the
    //! Shubnikov tags).  One hook, not "build, then fold": a structure whose build already searches every point's image under every op
    //! (a molecule's, to verify the mesh) folds from that search instead of repeating it (lane H 2026-11-25, SF₆'s mesh pilot).
    virtual FoldedMesh MakeFolded(const MeshParams& mp, const std::vector<Symmetry::SymOp>& ops) const=0;
};

} // namespace qchem::qcMesh
