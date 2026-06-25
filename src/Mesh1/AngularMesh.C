// File: AngularMesh.C  Abstract angular quadrature mesh + typed factory.
//
// Unit directions Omega_i and weights w_i normalised so that  sum_i w_i = 4*pi  (i.e.
// sum_i w_i g(Omega_i) ~ integral g dOmega).  Storage lives in the base; concretes fill it.
module;
#include <memory>
export module qchem.Mesh1.Angular;
export import qchem.Types;
export import qchem.Mesh1;            // AngularKind, MeshParams

export class AngularMesh
{
public:
    virtual ~AngularMesh() = default;
    const rvec3vec_t& Dirs() const {return itsD;}   //!< unit directions Omega_i
    const rvec_t&     W   () const {return itsW;}   //!< weights, sum = 4*pi
    size_t          size  () const {return itsD.size();}
protected:
    rvec3vec_t itsD;
    rvec_t     itsW;
};

//! \brief Build an angular mesh of the requested kind from the typed parameters.
export std::unique_ptr<AngularMesh> MakeAngular(const MeshParams&);
