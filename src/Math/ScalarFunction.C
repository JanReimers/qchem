// File: ScalarFunction.C  Mixin interface for real-space scalar functions.
//
// Pointwise only -- NO Mesh dependency.  (The old operator()(const Mesh&)/Gradient(const Mesh&)
// overloads were the ISP sin that dragged Mesh into every field; they were used only by the
// retired MeshIntegrator.  The free-function quadrature now lives in qcMesh.)
module;
export module qchem.ScalarFunction;
export import qchem.Types;

export namespace qchem
{
template <class T> class ScalarFunction
{
public:
    virtual ~ScalarFunction()  {};

    virtual T         operator()(const rvec3_t&) const=0;
    virtual vec3_t<T> Gradient  (const rvec3_t&) const=0;

    //! \brief Batch evaluation: \f$f(r)\f$ at every point at once.  The default just loops the pointwise
    //! \c op(); an overrider that can beat \f$O(N_\text{pts})\f$ pointwise cost (e.g. an FFT-backed field:
    //! inverse-transform once, then apply the functional) does so here.  The fitter samples a field on its
    //! own quadrature mesh via \c f(mesh.Points()) -- so the points are always that mesh, by construction.
    virtual vec_t<T> operator()(const rvec3vec_t& rs) const
    {
        vec_t<T> v(rs.size());
        for (size_t q=0;q<rs.size();q++) v[q]=(*this)(rs[q]);
        return v;
    }
};
}
