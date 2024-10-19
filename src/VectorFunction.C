// File: VectorFunction.C  Mixin interface for real space vector functions.



#include <VectorFunction.H>
#include "Mesh/Mesh.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include <cassert>
#include <complex>

template <class T> typename VectorFunction<T>::Mat VectorFunction<T>::operator() (const Mesh& mesh) const
{
    index_t n=GetVectorSize();
    Mat ret(MatLimits(n,mesh.GetNumPoints()));
    Fill(ret,T(0.0));
    Eval(mesh,ret);
    return ret;
}

template <class T> typename VectorFunction<T>::Vec3Mat VectorFunction<T>::Gradient(const Mesh& mesh) const
{
    index_t n=GetVectorSize();
    Vec3Mat ret(MatLimits(n,mesh.GetNumPoints()));
    Fill(ret,Vec3(0,0,0));
    EvalGrad(mesh,ret);
    return ret;
}

template <class T> void VectorFunction<T>::Eval(const Mesh& mesh, Mat& mat) const
{
    assert(mat.GetNumRows()==GetVectorSize    ());
    assert(mat.GetNumCols()==mesh.size());
    index_t n=GetVectorSize();
    int i=1;
    for (auto rw:mesh) 
    {
        const Vec& v((*this)(r(rw)));
        for (int j=1; j<=n; j++) mat(j,i) += v(j);
        i++;
    }
}

template <class T> void VectorFunction<T>::EvalGrad(const Mesh& mesh, Vec3Mat& mat) const
{
    assert(mat.GetNumRows()==GetVectorSize    ());
    assert(mat.GetNumCols()==mesh.size());
    index_t n=GetVectorSize();
    int i=1;
    for (auto rw:mesh) 
    {
        const Vec3Vec& v(Gradient(r(rw)));
        for (int j=1; j<=n; j++) mat(j,i) += v(j);
        i++;
    }
}

template class VectorFunction<double>;
template class VectorFunction<std::complex<double> >;
