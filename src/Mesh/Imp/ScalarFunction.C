// File: ScalarFunctionImp.C  Mixin interface for real space functions.
module;
#include <complex>

module qchem.ScalarFunction;
import qchem.Vector3D;


template <class T> vec_t<T> ScalarFunction<T>::operator() (const Mesh& mesh) const
{
    vec_t<T> v(mesh.size(),T(0.0));
    auto i(v.begin());
    for (auto rw:mesh) (*i++) += (*this)(r(rw));
    return v;
}

template <class T> vec3vec_t<T> ScalarFunction<T>::Gradient(const Mesh& mesh) const
{
    vec3vec_t<T> v(mesh.size(),rvec3_t(0,0,0));
    auto i(v.begin());
    for (auto rw:mesh) (*i++) += this->Gradient(r(rw));
    return v;
}

template class ScalarFunction<double>;
template class ScalarFunction<std::complex<double> >;