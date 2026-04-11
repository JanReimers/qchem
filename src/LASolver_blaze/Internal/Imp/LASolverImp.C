// File: LASolver_blaze/Internal/Imp/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigne solutions.
module;
#include <cassert>
#include <iostream>
#include <cmath>
#include "blaze/Math.h" 
module qchem.LASolver.Internal.Common;

using std::cout;
using std::endl;
//-------------------------------------------------------------------------
//
//  Common level.
//

template <class T> rsmat_t  LASolverCommon<T>::Transform(const rsmat_t& M) const
{
    mat_t<T> Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return MakeSymmetric(Mprime,"Test matrix");
}
template <class T> mat_t<T>  LASolverCommon<T>::Transform(const mat_t<T>& M) const
{
    mat_t<T> Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return Mprime;
}

template <class T> void LASolverCommon<T>::Rescale(mat_t<T>& V,const rvec_t& w)
{
    for (size_t j=0;j<V.columns();j++)
        column(V,j)/=sqrt(w[j]);
        
}

template <class T> void LASolverCommon<T>::Rescale(mat_t<T>& U,const rvec_t& s, mat_t<T>& Vt)
{
    for (size_t j=0;j<U.columns();j++)
        column(U,j)/=sqrt(s[j]);
    for (size_t i=0;i<Vt.rows();i++)
        row(Vt,i)/=sqrt(s[i]);
}

//
//  UsV from OML SVD which returns a diag matrix for s.
//
template <class T>  void LASolverCommon<T>::Truncate(mat_t<T>& U, rvec_t& s, mat_t<T>& Vt, double tol)
{
    assert(isSquare(U ));
    assert(isSquare(Vt));
    assert(U .rows()==s.size());
    assert(Vt.rows()==s.size());
    //  Find the index to truncate at.
    size_t index=0;
    for (auto i:s) 
        if (i>=tol)
            index++;
        else
            break;
    size_t n=s.size();
    assert(s[index-1]>=tol);
    assert(index==n || s[index]<tol);
    if (n-index >0)
        std::cout << "LASolverCommon truncating " << n-index << " singular values." << 
        " Min(s)="<< s[n-1] << " tol=" << tol << std::endl;
    //
    //  Two sets of vector limits
    //
    // VecLimits vl (U.GetRowLimits());
    // VecLimits vlt(vl.Low,index);
    // assert(vl.Low==1);
    //
    //  Truncate
    //
    if (index<n)
    {
        U =blaze::submatrix(U, 0,0,n    ,index);
        s =blaze::subvector(s ,0  ,index);
        Vt=blaze::submatrix(Vt,0,0,index,n    );

    }
    // U =U.SubMatrix(MatLimits(vl,vlt));
    // s =s.SubVector(vlt);
    // Vt=Vt.SubMatrix(MatLimits(vlt,vl));
}

//
//  Version for eigen routines which conventionally return ascending eigen vales.
//
template <class T>  void LASolverCommon<T>::Truncate(mat_t<T>& U, rvec_t& w,double tol)
{
    assert(U.columns()==w.size());
    //  Find the index to truncate at.
    size_t index=0;
    for (auto i:w) 
        if (i<tol)
            index++;
        else
            break;

    size_t n=w.size();
    assert(w[index]>=tol);
    if (index>0)
    {
        assert(w[index-1]<tol);
        std::cout << "LASolverCommon truncating " << index << " eigen values." << 
            " Min(w)="<< w[0] << " tol=" << tol << std::endl;
        //
        //  Two sets of vector limits
        //
        assert(U.rows()==U.columns());
        size_t nr=U.rows();
        // VecLimits vlt(index,nr);
        //
        //  Truncate
        // 
        U=blaze::submatrix( U, 0, index, nr, nr-index );
        // U =U .SubMatrix(MatLimits(vl,vlt));
        // U.ReBase(1,1);
        w=blaze::subvector(w,index,nr-index);
        // w =w .SubVector(vlt);
        // w.ReBase(1);
    }
}

template <class T>  smat_t<T> LASolverCommon<T>::MakeSymmetric(mat_t<T>& A,std::string name)
{
#ifdef false
    double del=::MakeSymmetric(A); // A=0.5*(A+~A)
    if (fabs(del) > 1e-9)
        std::cerr << "Warning: " << name << " asymmetry = " << del << " is big!" << std::endl;
#else
    A=0.5*(A+trans(A));
#endif
    return A;
}

 template class LASolverCommon<double>;



