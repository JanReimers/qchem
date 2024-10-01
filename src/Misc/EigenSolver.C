// File: EigenSolver.C  General eigen solver.

#include "Misc/EigenSolver.H"
#include "oml/diagonalmatrix.h"
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

//
//  UsV from Lapack SVD which returns a diag matrix for s.
//
template <class T>  void EigenSolverCommon<T>::Truncate(Mat& U, DMat& s, Mat& Vt, double tol)
{
    assert(U.GetLimits()==s.GetLimits());
    assert(U.GetLimits()==Vt.GetLimits());
    //  Find the index to truncate at.
    size_t index=0;
    for (auto i:s.GetDiagonal()) 
        if (i>=tol)
            index++;
        else
            break;
    assert(s(index,index)>=tol);
    assert(s(index+1,index+1)<tol);
    size_t n=s.GetDiagonal().size();
    std::cout << "EigenSolverCommon truncating " << n-index << " singular values." << 
    " Min(s)="<< s(n,n) << " yol=" << tol << std::endl;
    //
    //  Two sets of vector limits
    //
    VecLimits vl (U.GetRowLimits());
    VecLimits vlt(vl.Low,index);
    assert(vl.Low==1);
    //
    //  Truncate
    //
    U =U .SubMatrix(MatLimits(vl,vlt));
    s =s .SubMatrix(MatLimits(vlt,vlt));
    Vt=Vt.SubMatrix(MatLimits(vlt,vl));
}

//
//  UsV from OML SVD which returns a diag matrix for s.
//
template <class T>  void EigenSolverCommon<T>::Truncate(Mat& U, RVec& s, Mat& V, double tol)
{
    assert(U.GetRowLimits()==s.GetLimits());
    assert(U.GetColLimits()==s.GetLimits());
    assert(U.GetLimits()==V.GetLimits());
    //  Find the index to truncate at.
    size_t index=0;
    for (auto i:s) 
        if (i>=tol)
            index++;
        else
            break;
    assert(s(index)>=tol);
    assert(s(index+1)<tol);
    size_t n=s.size();
    std::cout << "EigenSolverCommon truncating " << n-index << " singular values." << 
    " Min(s)="<< s(n) << " yol=" << tol << std::endl;
    //
    //  Two sets of vector limits
    //
    VecLimits vl (U.GetRowLimits());
    VecLimits vlt(vl.Low,index);
    assert(vl.Low==1);
    //
    //  Truncate
    //
    U =U.SubMatrix(MatLimits(vl,vlt));
    s =s.SubVector(vlt);
    V =V.SubMatrix(MatLimits(vl,vlt));
}

//
//  Version for eigen routines which conventionally return ascending eigen vales.
//
template <class T>  void EigenSolverCommon<T>::Truncate(Mat& U, RVec& w,double tol)
{
    assert(U.GetColLimits()==w.GetLimits());
    //  Find the index to truncate at.
    size_t index=1;
    for (auto i:w) 
        if (i<tol)
            index++;
        else
            break;

    assert(w(index)>=tol);
    if (index>1)
    {
        assert(w(index-1)<tol);
        std::cout << "EigenSolverCommon truncating " << index-1 << " eigen values." << 
            " Min(w)="<< w(1) << " yol=" << tol << std::endl;
        //
        //  Two sets of vector limits
        //
        VecLimits vl (U.GetRowLimits());
        VecLimits vlt(index,vl.High);
        assert(vl.Low==1);
        //
        //  Truncate
        //
        U =U .SubMatrix(MatLimits(vl,vlt));
        U.ReBase(1,1);
        w =w .SubVector(vlt);
        w.ReBase(1);
    }
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

