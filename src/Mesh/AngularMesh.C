// File: AngularMesh.C  Angular quadrature mesh (concrete value type) + typed factory.
//
// Unit directions Omega_i and weights w_i normalised so that  sum_i w_i = 4*pi.  The schemes
// (Gauss / GaussLegendre / EulerMaclaren) are plain builder FUNCTIONS that return an AngularMesh --
// no class hierarchy.
module;
#include <utility>
export module qchem.Mesh.Angular;
export import qchem.Types;
export import qchem.Mesh;            // AngularKind, MeshParams

export namespace qchem::qcMesh
{

class AngularMesh
{
public:
    AngularMesh() = default;
    AngularMesh(rvec3vec_t d, rvec_t w) : itsD(std::move(d)), itsW(std::move(w)) {}
    const rvec3vec_t& Dirs() const {return itsD;}   //!< unit directions Omega_i
    const rvec_t&     W   () const {return itsW;}   //!< weights, sum = 4*pi
    size_t          size  () const {return itsD.size();}
private:
    rvec3vec_t itsD;
    rvec_t     itsW;
};

//! \brief Build an angular mesh of the requested kind from the typed parameters.
AngularMesh MakeAngular(const MeshParams&);

} //export namespace qchem::qcMesh

// Per-scheme builders -- declared NON-exported here, implemented in separate files; only MakeAngular
// (above) calls them.
namespace qchem::qcMesh
{
    AngularMesh GaussAngular        (int numDir);
    AngularMesh GaussLegendreAngular(int L);
    AngularMesh EulerMaclarenAngular(int L, int m);
}
