// File: Structure/MeshFactory.C  THE STRUCTURE MESH FACTORY: one MakeMesh overload per concrete structure (V-MESH-FACTORY, doc/MeshFactoryPlan.md §3).
//
// Which quadrature a geometry gets is decided HERE, by overload resolution on the concrete structure type: no runtime switch on a kind,
// no cast, and no library that is not the structure's ever switches on a recipe kind (the R2.17 ruling, UnitCell.C history: "a BASIS was
// switching on cellKind to decide how a STRUCTURE builds its own mesh").  The bodies are the geometry's own builders, kept in qcStructure
// beside the geometry they read: the atom-centred Becke grid (Atom, Molecule), the periodic Becke grid and the uniform raster (UnitCell),
// and the site-adapted invariant build under an imposed op set.
//
// A NEW STRUCTURE KIND (a slab) adds an overload here, in a non-hub module (OCP: adding, not modifying).  A library that must choose a
// recipe itself holds a qcMesh::MeshMaker (qchem.Mesh.Maker); StructureMeshMaker<S> below is its one realization, built where S is known.
//
// TRANSITIONAL (§8 step 2): Structure::CreateIntegrationMesh and its overrides still exist and DELEGATE here; §8 step 4 deletes them.
module;
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
export module qchem.Structure.MeshFactory;
export import qchem.Mesh.Maker;        // qcMesh::MeshMaker, Mesh, MeshParams
export import qchem.Structure;         // Atom, Molecule
export import qchem.UnitCell;          // UnitCell
export import qchem.Symmetry.SymOp;    // Symmetry::SymOp (the imposed overload's op set)

export namespace qchem
{

//! A single atom's mesh: its atom-centred radial x angular product grid (the Becke weight is identically 1 with no other centre).
qcMesh::Mesh MakeMesh(const Atom&,     const qcMesh::MeshParams&);
//! A molecule's mesh: the multi-centre Becke-partitioned grid, one site block per atom (Becke 1988).
qcMesh::Mesh MakeMesh(const Molecule&, const qcMesh::MeshParams&);
//! A periodic cell's mesh; \c mp.cellKind selects it: \c Uniform, the midpoint raster (\c mp.eCut or \c mp.nUniform per axis); \c Becke,
//! the atom-centred periodic fuzzy-Voronoi quadrature (image series converged by magnitude screening, no radius anywhere: pin 1).
qcMesh::Mesh MakeMesh(const UnitCell&, const qcMesh::MeshParams&);
//! The same cell mesh made INVARIANT under \a ops (§6a imposed symmetry): site-adapted by construction for \c Becke (W2b), the generic
//! raster group-averaged onto the op orbits for \c Uniform.  \a ops non-empty.
qcMesh::Mesh MakeMesh(const UnitCell&, const qcMesh::MeshParams&, const std::vector<Symmetry::SymOp>& ops);
//! An atom's mesh invariant under \a ops: the radial-only set (\c AngularKind::None), whose angular part is analytic, so ANY rotation about
//! the atom leaves it invariant.  An angular recipe THROWS: its finite direction set would claim a symmetry it does not carry.
qcMesh::Mesh MakeMesh(const Atom&,     const qcMesh::MeshParams&, const std::vector<Symmetry::SymOp>& ops);
//! A molecule's mesh invariant under a point group: NOT BUILT YET (OOD-SOLID V-FOLD-NEUTRAL) -- THROWS by name.  A molecule today is only
//! ever handed \c qcMesh::NoSymmetry, whose empty op set never reaches this.
qcMesh::Mesh MakeMesh(const Molecule&, const qcMesh::MeshParams&, const std::vector<Symmetry::SymOp>& ops);

//! The ORBIT FOLD of a cell mesh invariant under \a ops: the orbit partition (fractional, torus-metric point matching), and on a
//! Shubnikov op set the per-op spin actions and the flip-fixed points (m = 0 there exactly).
qcMesh::FoldedMesh FoldOrbits(const UnitCell&, qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>& ops);
//! An atom's radial-only mesh: one direction per shell, so every orbit is a single point -- the trivial fold.
qcMesh::FoldedMesh FoldOrbits(const Atom&,     qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>& ops);
//! A molecule's: NOT BUILT YET (OOD-SOLID V-FOLD-NEUTRAL) -- THROWS by name.
qcMesh::FoldedMesh FoldOrbits(const Molecule&, qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>& ops);

//! \brief THE concrete \c qcMesh::MeshMaker: closes over one structure of concrete type \a S and builds through \c MakeMesh(const S&, mp).
//! One class template instead of a hand-written hierarchy parallel to the structures (§2 (ii)).  The structure is DEMANDED: no null.
template <class S> class StructureMeshMaker : public qcMesh::MeshMaker
{
public:
    explicit StructureMeshMaker(std::shared_ptr<const S> s) : itsStructure(std::move(s))
    {
        if (!itsStructure) throw std::invalid_argument("StructureMeshMaker: a null structure -- a mesh maker needs the geometry it builds over");
    }
    using qcMesh::MeshMaker::Make;   // the symmetric overload, the face's template method
    qcMesh::Mesh Make(const qcMesh::MeshParams& mp) const override {return MakeMesh(*itsStructure, mp);}
protected:
    qcMesh::Mesh MakeInvariant(const qcMesh::MeshParams& mp, const std::vector<Symmetry::SymOp>& ops) const override
        {return MakeMesh(*itsStructure, mp, ops);}
    qcMesh::FoldedMesh FoldOrbits(qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>& ops) const override
        {return qchem::FoldOrbits(*itsStructure, std::move(mesh), ops);}
private:
    std::shared_ptr<const S> itsStructure;
};

//! A maker over a COPY of \a s (a concrete \c Atom, \c Molecule or \c UnitCell): for a caller that holds its geometry by value -- a
//! facade's point cell, a unit test's -- and hands the maker to a library that builds at a recipe only it knows (the GPW basis, §4.4).
template <class S> std::shared_ptr<const qcMesh::MeshMaker> MakeMeshMaker(const S& s)
{
    return std::make_shared<const StructureMeshMaker<S>>(std::make_shared<const S>(s));
}

} // namespace qchem
