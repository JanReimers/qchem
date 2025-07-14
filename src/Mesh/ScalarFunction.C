// File: ScalarFunction.C  Mixin interface for real space functions.



#include "oml/vector.h"
#include <cassert>
#include <iostream>
#include <complex>
#include <valarray>
#include <Mesh/ScalarFunction.H>


template <class T> typename ScalarFunction<T>::Vec ScalarFunction<T>::operator() (const Mesh& mesh) const
{
    Vec v(mesh.size());
    Fill(v,T(0.0));
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

template <class T> typename ScalarFunction<T>::Vec3Vec ScalarFunction<T>::Gradient(const Mesh& mesh) const
{
    Vec3Vec v(mesh.size());
    Fill(v,Vec3(0,0,0));
    auto i(v.begin());
    for (auto rw:mesh) (*i++) += this->Gradient(r(rw));
    return v;
}

template class ScalarFunction<double>;
template class ScalarFunction<std::complex<double> >;
