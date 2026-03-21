// File: LASolver_blaze/Internal/Imp/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigne solutions.
module;
#include <cassert>
#include <iostream>
#include <cmath>
#include "blaze/Math.h" 
module qchem.LASolver_blaze.Internal.Common;
//-------------------------------------------------------------------------
//
//  Common level.
//

template <class T> typename LASolver_blaze<T>::rsmat_t  LASolverCommon_blaze<T>::Transform(const rsmat_t& M) const
{
    mat_t Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return MakeSymmetric(Mprime,"Test matrix");
}
template <class T> typename LASolver_blaze<T>::mat_t  LASolverCommon_blaze<T>::Transform(const mat_t& M) const
{
    mat_t Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return Mprime;
}

template <class T> void LASolverCommon_blaze<T>::Rescale(mat_t& V,const rvec_t& w)
{
    for (size_t j=0;j<V.columns();j++)
        column(V,j)/=sqrt(w[j]);
        
}

template <class T> void LASolverCommon_blaze<T>::Rescale(mat_t& U,const rvec_t& s, mat_t& Vt)
{
    for (size_t j=0;j<U.columns();j++)
        column(U,j)/=sqrt(s[j]);
    for (size_t i=0;i<Vt.rows();i++)
        row(Vt,i)/=sqrt(s[i]);
}

//
//  UsV from OML SVD which returns a diag matrix for s.
//
template <class T>  void LASolverCommon_blaze<T>::Truncate(mat_t& U, rvec_t& s, mat_t& Vt, double tol)
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
    assert(s[index]>=tol);
    assert(index==n || s[index+1]<tol);
    if (n-index >0)
        std::cout << "LASolverCommon truncating " << n-index << " singular values." << 
        " Min(s)="<< s[n] << " tol=" << tol << std::endl;
    //
    //  Two sets of vector limits
    //
    // VecLimits vl (U.GetRowLimits());
    // VecLimits vlt(vl.Low,index);
    // assert(vl.Low==1);
    //
    //  Truncate
    //
    U =blaze::submatrix(U, 0,n,0,index);
    s =blaze::subvector(s ,0,index);
    Vt=blaze::submatrix(Vt,0,index,0,n);
    // U =U.SubMatrix(MatLimits(vl,vlt));
    // s =s.SubVector(vlt);
    // Vt=Vt.SubMatrix(MatLimits(vlt,vl));
}

//
//  Version for eigen routines which conventionally return ascending eigen vales.
//
template <class T>  void LASolverCommon_blaze<T>::Truncate(mat_t& U, rvec_t& w,double tol)
{
    assert(U.columns()==w.size());
    //  Find the index to truncate at.
    size_t index=0;
    for (auto i:w) 
        if (i<tol)
            index++;
        else
            break;

    assert(w[index]>=tol);
    if (index>0)
    {
        assert(w[index]<tol);
        std::cout << "LASolverCommon truncating " << index+1 << " eigen values." << 
            " Min(w)="<< w[0] << " yol=" << tol << std::endl;
        //
        //  Two sets of vector limits
        //
        assert(U.rows()==U.columns());
        size_t nr=U.rows();
        // VecLimits vlt(index,nr);
        //
        //  Truncate
        //
        U=blaze::submatrix( U, 0UL, nr, index, nr );
        // U =U .SubMatrix(MatLimits(vl,vlt));
        // U.ReBase(1,1);
        w=blaze::subvector(w,index,nr);
        // w =w .SubVector(vlt);
        // w.ReBase(1);
    }
}

template <class T>  LASolver_blaze<T>::smat_t LASolverCommon_blaze<T>::MakeSymmetric(mat_t& A,std::string name)
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

 template class LASolverCommon_blaze<double>;



