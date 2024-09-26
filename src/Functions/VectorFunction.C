// File: VectorFunction.C  Mixin interface for real space vector functions.



#include "Functions/VectorFunction.H"
#include "Mesh/Mesh.H"
#include "Mesh/MeshBrowser.H"
#include "oml/matrix.h"
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
    assert(mat.GetNumCols()==mesh.GetNumPoints());
    index_t n=GetVectorSize();
    MeshBrowser b(mesh);
    typename Mat::Subscriptor s(mat);
    for (int i=1; b; b++,i++)
    {
        const Vec& v(operator()(b.R()));
        for (int j=1; j<=n; j++) s(j,i) += v(j);
    }
}

template <class T> void VectorFunction<T>::EvalGrad(const Mesh& mesh, Vec3Mat& mat) const
{
    assert(mat.GetNumRows()==GetVectorSize    ());
    assert(mat.GetNumCols()==mesh.GetNumPoints());
    index_t n=GetVectorSize();
    MeshBrowser b(mesh);
    typename Vec3Mat::Subscriptor s(mat);
    for (int i=1; b; b++,i++)
    {
        const Vec3Vec& v(Gradient(b.R()));
        for (int j=1; j<=n; j++) s(j,i) += v(j);
    }
}

template class VectorFunction<double>;
template class VectorFunction<std::complex<double> >;
