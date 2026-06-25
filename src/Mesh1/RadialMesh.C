// File: RadialMesh.C  Abstract radial quadrature mesh + typed factory.
//
// A 1D radial mesh: nodes r_i and weights w_i with the r^2 jacobian ALREADY folded into w_i, so
// that  sum_i w_i f(r_i) ~ integral_0^inf r^2 f(r) dr.  Storage lives in the base (every concrete
// stores the same two arrays); the concretes only differ in how they fill them.
module;
#include <memory>
export module qchem.Mesh1.Radial;
export import qchem.Types;
export import qchem.Mesh1;            // RadialKind, MeshParams

export namespace qcMesh1
{

class RadialMesh
{
public:
    virtual ~RadialMesh() = default;
    const rvec_t& R() const {return itsR;}   //!< r_i
    const rvec_t& W() const {return itsW;}   //!< w_i  (includes the r^2 jacobian)
    size_t      size() const {return itsR.size();}
protected:
    rvec_t itsR, itsW;
};

//! \brief Build a radial mesh of the requested kind from the typed parameters.
std::unique_ptr<RadialMesh> MakeRadial(const MeshParams&);

} //export namespace qcMesh1
