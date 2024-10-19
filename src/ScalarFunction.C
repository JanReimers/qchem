// File: ScalarFunction.C  Mixin interface for real space functions.



#include <ScalarFunction.H>
#include "Mesh/Mesh.H"
#include "Mesh/MeshBrowser.H"
#include "oml/vector.h"
#include <cassert>
#include <iostream>
#include <complex>

template <class T> typename ScalarFunction<T>::Vec ScalarFunction<T>::operator() (const Mesh& mesh) const
{
    Vec ret(mesh.GetNumPoints());
    Fill(ret,T(0.0));
    Eval(mesh, ret);
    return ret;
}

template <class T> typename ScalarFunction<T>::Vec3Vec ScalarFunction<T>::Gradient(const Mesh& mesh) const
{
    Vec3Vec ret(mesh.GetNumPoints());
    Fill(ret,Vec3(0,0,0));
    EvalGrad(mesh, ret);
    return ret;
}

//
//  The eval funciont assume that v is already initialized.  New results are just added on.
//
template <class T> void ScalarFunction<T>::Eval(const Mesh& mesh, Vec& v) const
{
    assert(v.size()==mesh.GetNumPoints());
    MeshBrowser          b(mesh);
    typename Vec::iterator i(v.begin());
    for (; b&&i!=v.end(); b++,i++)
        *i+=this->operator()(b.R());
}

template <class T> void ScalarFunction<T>::EvalGrad(const Mesh& mesh, Vec3Vec& v) const
{
    assert(v.size()==mesh.GetNumPoints());
    MeshBrowser            b(mesh);
    typename Vec3Vec::iterator i(v.begin());
    for (; b&&i!=v.end(); b++,i++) *i+=this->Gradient(b.R());
}

template class ScalarFunction<double>;
template class ScalarFunction<std::complex<double> >;
