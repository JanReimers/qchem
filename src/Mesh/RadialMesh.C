// File: RadialMesh.C  Radial quadrature mesh (concrete value type) + typed factory.
//
// A 1D radial mesh: nodes r_i and weights w_i with the r^2 jacobian ALREADY folded into w_i, so
// that  sum_i w_i f(r_i) ~ integral_0^inf r^2 f(r) dr.  The schemes (MHL/Log/Linear) are plain
// builder FUNCTIONS that return a RadialMesh -- no class hierarchy (the old per-scheme classes did
// nothing but fill these two arrays).
module;
#include <utility>
export module qchem.Mesh.Radial;
export import qchem.Types;
export import qchem.Mesh;            // RadialKind, MeshParams

export namespace qcMesh
{

class RadialMesh
{
public:
    RadialMesh() = default;
    RadialMesh(rvec_t r, rvec_t w) : itsR(std::move(r)), itsW(std::move(w)) {}
    const rvec_t& R() const {return itsR;}   //!< r_i
    const rvec_t& W() const {return itsW;}   //!< w_i  (includes the r^2 jacobian)
    size_t      size() const {return itsR.size();}
private:
    rvec_t itsR, itsW;
};

//! \brief Build a radial mesh of the requested kind from the typed parameters.
RadialMesh MakeRadial(const MeshParams&);

} //export namespace qcMesh

// Per-scheme builders -- declared NON-exported here, implemented in separate files; only MakeRadial
// (above) calls them.
namespace qcMesh
{
    RadialMesh MHLRadial   (int NumPoints, int m, double alpha);
    RadialMesh LogRadial   (double start, double stop, int NumPoints);
    RadialMesh LinearRadial(double start, double stop, int NumPoints);
}
