// File: ScalarFunctionImp.C  Mixin interface for real space functions.
module;
#include <cassert>
#include <iostream>
#include <complex>
#include <valarray>

module qchem.ScalarFunction;
import oml.Vector3D;


template <class T> vec_t<T> ScalarFunction<T>::operator() (const Mesh& mesh) const
{
    vec_t<T> v(mesh.size(),T(0.0));
    auto i(v.begin());
    for (auto rw:mesh) (*i++) += (*this)(r(rw));
    return v;
}

template <class T> typename ScalarFunction<T>::va_t ScalarFunction<T>::operator() (const rva_t& rs,RVec3 dir) const
{
    assert(norm(dir)==1.0);
    va_t v(T(0),rs.size());
    auto i(std::begin(v));
    for (auto r:rs) (*i++) += (*this)(r*dir);
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