// File: VectorFunction.C  Mixin interface for real-space vector (basis) functions.
//
// Pointwise only -- NO Mesh dependency.  When evaluated at a point r it returns a vector of values
// [phi_i(r)].  (The old operator()(const Mesh&)/Gradient(const Mesh&) overloads -- the ISP sin --
// were used only by the retired MeshIntegrator; quadrature now lives in qcMesh.)
module;
export module qchem.VectorFunction;
export import qchem.Types;

namespace qchem {

export template <class T> class VectorFunction
{
public:
    virtual ~VectorFunction()  {};

    virtual size_t       GetVectorSize()             const=0;
    virtual vec_t<T>     operator() (const rvec3_t&) const=0;
    virtual vec3vec_t<T> Gradient   (const rvec3_t&) const=0;
};

} // namespace qchem