// File: LASolver.C  General eigen solver.

#include "Imp/LASolver/LASolverImp.H"
#include "oml/diagonalmatrix.h"
#include "oml/vector.h"
#include <cmath>

//-------------------------------------------------------------------------
//
//  Common level.
//
template <class T> void LASolverCommon<T>::Rescale(Mat& V,const RVec& w)
{
    for (auto j:V.cols())
        V.GetColumn(j)/=sqrt(w(j));
        
}

template <class T> void LASolverCommon<T>::Rescale(Mat& U,const RVec& s, Mat& Vt)
{
    for (auto j:U.cols())
        U.GetColumn(j)/=sqrt(s(j));
    for (auto i:Vt.rows())
        Vt.GetRow(i)/=sqrt(s(i));
        
}

//
//  UsV from OML SVD which returns a diag matrix for s.
//
template <class T>  void LASolverCommon<T>::Truncate(Mat& U, RVec& s, Mat& Vt, double tol)
{
    assert(U.GetRowLimits()==s.GetLimits());
    assert(U.GetColLimits()==s.GetLimits());
    assert(U.GetLimits()==Vt.GetLimits());
    //  Find the index to truncate at.
    size_t index=0;
    for (auto i:s) 
        if (i>=tol)
            index++;
        else
            break;
    size_t n=s.size();
    assert(s(index)>=tol);
    assert(index==n || s(index+1)<tol);
    std::cout << "LASolverCommon truncating " << n-index << " singular values." << 
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
    Vt=Vt.SubMatrix(MatLimits(vlt,vl));
}

//
//  Version for eigen routines which conventionally return ascending eigen vales.
//
template <class T>  void LASolverCommon<T>::Truncate(Mat& U, RVec& w,double tol)
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
        std::cout << "LASolverCommon truncating " << index-1 << " eigen values." << 
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

template class LASolver<double>;
template class LASolverCommon<double>;


#include "Imp/LASolver/LASolverOML.H"
#include "Imp/LASolver/LASolverLapack.H"


template <class T> LASolver<T>* LASolver<T>::
    Factory(const LinearAlgebraParams& lap)
{
    LASolver<T>* ret=0;
    switch (lap.LinearAlgebraPackage)
    {
    case qchem::OML:
        switch (lap.BasisOrthoAlgorithm)
        {
        case qchem::Cholsky :
            ret=new LASolverOMLCholsky<T>(lap);
            break;
        case qchem::Eigen :
            ret=new LASolverOMLEigen<T>(lap);
            break;
        case qchem::SVD :
            ret=new LASolverOMLSVD<T>(lap);
            break;
        }
        break;
   case qchem::Lapack:
        switch (lap.BasisOrthoAlgorithm)
        {
        case qchem::Cholsky :
            ret=new LASolverLapackCholsky<T>(lap);
            break;
        case qchem::Eigen :
            ret=new LASolverLapackEigen<T>(lap);
            break;
        case qchem::SVD :
            ret=new LASolverLapackSVD<T>(lap);
            break;
        break;
        }
    }
    assert(ret);
    return ret;
}

