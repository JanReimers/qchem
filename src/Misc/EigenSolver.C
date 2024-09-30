// File: EigenSolver.C  General eigen solver.

#include "Misc/EigenSolver.H"
#include "oml/vector.h"
#include <cmath>

//-------------------------------------------------------------------------
//
//  Common level.
//
template <class T> void EigenSolverCommon<T>::Rescale(Mat& V,const RVec& w)
{
    for (auto j:V.cols())
        V.GetColumn(j)/=sqrt(w(j));
        
}

template <class T> void EigenSolverCommon<T>::Rescale(Mat& U,const RVec& s, Mat& Vt)
{
    for (auto j:U.cols())
        U.GetColumn(j)/=sqrt(s(j));
    for (auto i:Vt.rows())
        Vt.GetRow(i)/=sqrt(s(i));
        
}


template class EigenSolver<double>;
template class EigenSolverCommon<double>;


#include "Misc/EigenSolverOML.H"
#include "Misc/EigenSolverLapack.H"


template <class T> EigenSolver<T>* EigenSolver<T>::
    Factory(EigenSolver<T>::Pkg pkg,EigenSolver<T>::Ortho ortho,const SMat& S, double tolerance)
{
    EigenSolver<T>* ret=0;
    switch (pkg)
    {
    case OML:
        switch (ortho)
        {
        case Cholsky :
            ret=new EigenSolverOMLCholsky<T>(S,tolerance);
            break;
        case Eigen :
            ret=new EigenSolverOMLEigen<T>(S,tolerance);
            break;
        case SVD :
            ret=new EigenSolverOMLSVD<T>(S,tolerance);
            break;
        }
        break;
   case Lapack:
        switch (ortho)
        {
        case Cholsky :
            ret=new EigenSolverLapackCholsky<T>(S,tolerance);
            break;
        case Eigen :
            ret=new EigenSolverLapackEigen<T>(S,tolerance);
            break;
        case SVD :
            ret=new EigenSolverLapackSVD<T>(S,tolerance);
            break;
        break;
        }
    }
    assert(ret);
    return ret;
}

