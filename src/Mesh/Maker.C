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
export module qchem.Mesh.Maker;
export import qchem.Mesh;          // qcMesh::Mesh, the built object Make returns
export import qchem.Mesh.Recipe;   // qcMesh::MeshParams, the recipe Make reads

export namespace qchem::qcMesh
{

//! \brief Builds a real-space integration mesh from a recipe over ONE fixed structure (the one the realization closes over).
class MeshMaker
{
public:
    virtual ~MeshMaker()=default;
    //! The mesh of the bound structure at recipe \a mp.  Throws what the structure's builder throws (a pin-30 refusal, a broken invariant).
    virtual Mesh Make(const MeshParams& mp) const=0;
};

} // namespace qchem::qcMesh
