// File: Imp/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigne solutions.
module;
#include <cassert>
#include <iostream>
#include <cmath>
module qchem.LASolver.Internal.Common;
//-------------------------------------------------------------------------
//
//  Common level.
//

template <class T> typename LASolver<T>::RSMat  LASolverCommon<T>::Transform(const RSMat& M) const
{
    Mat Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return MakeSymmetric(Mprime,"Test matrix");
}
template <class T> typename LASolver<T>::Mat  LASolverCommon<T>::Transform(const Mat& M) const
{
    Mat Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return Mprime;
}

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
    if (n-index >0)
        std::cout << "LASolverCommon truncating " << n-index << " singular values." << 
        " Min(s)="<< s(n) << " tol=" << tol << std::endl;
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

template <class T>  typename LASolverCommon<T>::SMat LASolverCommon<T>::MakeSymmetric(Mat& A,std::string name)
{
#ifdef false
    double del=::MakeSymmetric(A); // A=0.5*(A+~A)
    if (fabs(del) > 1e-9)
        std::cerr << "Warning: " << name << " asymmetry = " << del << " is big!" << std::endl;
#else
    ::MakeSymmetric(A); // A=0.5*(A+~A)
#endif
    return A;
}

 template class LASolverCommon<double>;



