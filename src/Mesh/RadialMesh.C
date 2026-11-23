// File: RadialMesh.C  The radial quadrature mesh's typed factory (the value type RadialMesh is in qchem.Mesh).
//
// A 1D radial mesh: nodes r_i and weights w_i with the r^2 jacobian ALREADY folded into w_i, so
// that  sum_i w_i f(r_i) ~ integral_0^inf r^2 f(r) dr.  The schemes (MHL/Log/Linear/DE2) are plain
// builder FUNCTIONS that return a RadialMesh -- no class hierarchy (the old per-scheme classes did
// nothing but fill these two arrays).
module;
#include <utility>
export module qchem.Mesh.Radial;
export import qchem.Types;
export import qchem.Mesh;            // RadialMesh (the built value type lives there, V-MESH-FACTORY step 5)
import qchem.Mesh.Recipe;            // MeshParams, RadialKind: this module is the recipe-side FACTORY of a RadialMesh

export namespace qchem::qcMesh
{


//! \brief Build a radial mesh of the requested kind from the typed parameters.
RadialMesh MakeRadial(const MeshParams&);

} //export namespace qchem::qcMesh

// Per-scheme builders -- declared NON-exported here, implemented in separate files; only MakeRadial
// (above) calls them.
namespace qchem::qcMesh
{
    RadialMesh MHLRadial   (int NumPoints, int m, double alpha);
    RadialMesh LogRadial   (double start, double stop, int NumPoints);
    RadialMesh LinearRadial(double start, double stop, int NumPoints);
    RadialMesh DE2Radial   (int NumPoints, double alpha, double rMin, double rMax);
}
