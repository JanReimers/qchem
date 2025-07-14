// File: VectorFunction.C  Mixin interface for real space vector functions.



#include "oml/matrix.h"
#include "oml/vector.h"
#include <cassert>
#include <complex>
#include <Mesh/VectorFunction.H>

template <class T> typename VectorFunction<T>::Mat VectorFunction<T>::operator() (const Mesh& mesh) const
{
    index_t n=GetVectorSize();
    Mat m(MatLimits(n,mesh.size()));
    Fill(m,T(0.0));
    int i=1;
    for (auto rw:mesh) 
    {
        const Vec& v((*this)(r(rw)));
        for (int j=1; j<=n; j++) m(j,i) += v(j);
        i++;
    }
    return m;
}

template <class T> typename VectorFunction<T>::Vec3Mat VectorFunction<T>::Gradient(const Mesh& mesh) const
{
    index_t n=GetVectorSize();
    Vec3Mat m(MatLimits(n,mesh.size()));
    Fill(m,Vec3(0,0,0));
    int i=1;
    for (auto rw:mesh) 
    {
        const Vec3Vec& v(Gradient(r(rw)));
        for (int j=1; j<=n; j++) m(j,i) += v(j);
        i++;
    }
    return m;
}

template class VectorFunction<double>;
template class VectorFunction<std::complex<double> >;
