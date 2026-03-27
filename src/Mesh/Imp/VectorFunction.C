// File: ScalarFunctionImp.C  Mixin interface for real space vector functions.
module;
#include <cassert>
#include <complex>
module qchem.VectorFunction;
import oml;

template <class T> typename VectorFunction<T>::Mat VectorFunction<T>::operator() (const Mesh& mesh) const
{
    size_t n=GetVectorSize();
    Mat m(n,mesh.size(),T(0.0));
    int i=0;
    for (auto rw:mesh) 
    {
        const Vec& v((*this)(r(rw)));
        for (int j=0; j<n; j++) m(j,i) += v[j];
        i++;
    }
    return m;
}

template <class T> typename VectorFunction<T>::Vec3Mat VectorFunction<T>::Gradient(const Mesh& mesh) const
{
    size_t n=GetVectorSize();
    Vec3Mat m(n,mesh.size(),Vec3(0,0,0));
    int i=0;
    for (auto rw:mesh) 
    {
        const Vec3Vec& v(Gradient(r(rw)));
        for (int j=0; j<n; j++) m(j,i) += v[j];
        i++;
    }
    return m;
}

template class VectorFunction<double>;
template class VectorFunction<std::complex<double> >;

