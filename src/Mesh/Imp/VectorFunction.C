// File: ScalarFunctionImp.C  Mixin interface for real space vector functions.
module;
#include <cassert>
#include <complex>
module qchem.VectorFunction;
import oml.Vector3D;

template <class T> mat_t<T> VectorFunction<T>::operator() (const Mesh& mesh) const
{
    size_t n=GetVectorSize();
    mat_t<T> m(n,mesh.size(),T(0.0));
    int i=0;
    for (auto rw:mesh) 
    {
        const vec_t<T>& v((*this)(r(rw)));
        for (int j=0; j<n; j++) m(j,i) += v[j];
        i++;
    }
    return m;
}

template <class T> vec3mat_t<T> VectorFunction<T>::Gradient(const Mesh& mesh) const
{
    size_t n=GetVectorSize();
    vec3mat_t<T> m(n,mesh.size(),vec3_t<T>(0,0,0));
    int i=0;
    for (auto rw:mesh) 
    {
        const vec3vec_t<T>& v(Gradient(r(rw)));
        for (int j=0; j<n; j++) m(j,i) += v[j];
        i++;
    }
    return m;
}

template class VectorFunction<double>;
template class VectorFunction<std::complex<double> >;

