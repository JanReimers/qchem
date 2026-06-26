// File: ScalarFunction.C  Mixin interface for real-space scalar functions.
//
// Pointwise only -- NO Mesh dependency.  (The old operator()(const Mesh&)/Gradient(const Mesh&)
// overloads were the ISP sin that dragged Mesh into every field; they were used only by the
// retired MeshIntegrator.  The free-function quadrature now lives in qcMesh1.)
module;
export module qchem.ScalarFunction;
export import qchem.Types;

export template <class T> class ScalarFunction
{
public:
    virtual ~ScalarFunction()  {};

    virtual T         operator()(const rvec3_t&) const=0;
    virtual vec3_t<T> Gradient  (const rvec3_t&) const=0;
};
