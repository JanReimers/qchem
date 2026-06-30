// File: Builder.C  Incremental Mesh accumulation, the efficient way.
//
// When a builder appends points one at a time and the final count is not known up front (the Becke
// molecular mesh: each atom contributes a variable number of surviving points), Mesh::Append's
// resize-by-one is O(N^2).  MeshBuilder accumulates into the project's doubling blazem::VecBuilder
// (O(N) total) and hands over a finished Mesh via the from-arrays constructor.  Opt-in module so
// the core qchem.Mesh stays free of the qchem.Blaze dependency.
module;
export module qchem.Mesh.Builder;
export import qchem.Mesh;
import qchem.Blaze;          // blazem::VecBuilder

export namespace qchem::qcMesh
{

class MeshBuilder
{
public:
    void   Append(const rvec3_t& r, double w) {itsR.Append(r); itsW.Append(w);}
    size_t size  () const                     {return itsW.size();}
    Mesh   take  ()                           {return Mesh(itsR.take(), itsW.take());}
private:
    blazem::VecBuilder<rvec3_t> itsR;
    blazem::VecBuilder<double>  itsW;
};

} //export namespace qchem::qcMesh
