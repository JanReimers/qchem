// File: ScalarFunction.C  Mixin interface for real-space scalar functions.
module;
export module qchem.ScalarFunction;
export import qchem.Types;
import qchem.PointSampling;   // PointSampling (SampleOn's argument)

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
    //! \brief Sample me on a POINT SET -- the double dispatch that lets the point set choose HOW (qchem.PointSampling).
    //! The default evaluates me at its \c Points(), bit-identical to \c f(points); a field held as a Fourier series
    //! overrides it to ask the point set for its sampler (a NUFFT on any set, an exact inverse FFT on the uniform
    //! raster), so neither side ever asks what kind the other is.
    virtual vec_t<T> SampleOn(const PointSampling& ps) const {return (*this)(ps.Points());}
};
}
