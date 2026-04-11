// File: LASolver/Internal/Imp/LASolverLapack.C  General eigen solver.
module;
#include <cassert>
#include <iostream>
#include "blaze/Math.h" 
module qchem.LASolver.Internal.Lapack;


template <class T> void LASolverEigen<T>::SetBasisOverlap(const smat_t<T>& S)
{
    rvec_t d;
    mat_t<T>  U;
    blaze::eigen(S,d,U);
    Truncate(U,d,itsTruncationTolerance);
    Rescale(U,d);
    LASolverCommon<T>::AssignVs(U,trans(U));
    LASolverCommon<T>::Diag=d; //Preserve eigend values.
}

template <class T> rsmat_t LASolverEigen<T>::Inverse(const rsmat_t& S) const
{
    rvec_t d;
    mat_t<T>  U;
    blaze::eigen(S,d,U);
    // auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    Truncate(U,d,itsTruncationTolerance);
    blaze::DiagonalMatrix<rmat_t> winv(d.size());
    blaze::diagonal(winv)=1.0/d;
    mat_t<T> Sfull(U*winv*trans(U));
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverEigen<T>::Rescale(mat_t<T>& V,const rvec_t& w)
{
    for (size_t j=0;j<V.columns();j++)
        column(V,j)/=sqrt(w[j]);
        
}


//
//  Version for eigen routines which conventionally return ascending eigen vales.
//
template <class T>  void LASolverEigen<T>::Truncate(mat_t<T>& U, rvec_t& w,double tol)
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


template <class T> void LASolverSVD<T>::SetBasisOverlap(const smat_t<T>& S)
{
    rvec_t s;
    mat_t<T>  U,Vt;
    blaze::svd(S,U,s,Vt);
    Truncate(U,s,Vt,itsTruncationTolerance);
    Rescale(U,s,Vt);
    LASolverCommon<T>::AssignVs(U,Vt);
    LASolverCommon<T>::Diag=s; //Preserve SVs.
}

template <class T> rsmat_t LASolverSVD<T>::Inverse(const rsmat_t& S) const
{
    rvec_t s;
    mat_t<T>  U,Vt;
    blaze::svd(S,U,s,Vt);
    Truncate(U,s,Vt,itsTruncationTolerance);
    blaze::DiagonalMatrix<rmat_t> sinv(s.size());
    blaze::diagonal(sinv)=1.0/s;
    mat_t<T>  Sfull(trans(Vt)*sinv*trans(U));
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverSVD<T>::Rescale(mat_t<T>& U,const rvec_t& s, mat_t<T>& Vt)
{
    for (size_t j=0;j<U.columns();j++)
        column(U,j)/=sqrt(s[j]);
    for (size_t i=0;i<Vt.rows();i++)
        row(Vt,i)/=sqrt(s[i]);
}

//
//  UsV from OML SVD which returns a diag matrix for s.
//
template <class T>  void LASolverSVD<T>::Truncate(mat_t<T>& U, rvec_t& s, mat_t<T>& Vt, double tol)
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

template <class T> void LASolverCholsky<T>::SetBasisOverlap(const smat_t<T>& S)
{
    mat_t<T> Sm(S);
    blaze::potrf( Sm, 'U' );
    blaze::trtri( Sm, 'U', 'N' );
    size_t N=Sm.rows();
    for (size_t i=0;i<N;i++)
        for (size_t j=i+1;j<N;j++)
            Sm(j,i)=0;
    blaze::UpperMatrix< mat_t<T>> V(Sm);
    

    LASolverCommon<T>::AssignVs(V,blaze::trans(V)); //V=U^-1, Vd=transpose(U^-1)
    LASolverCommon<T>::Diag=blaze::diagonal(V);
}

template <class T> rsmat_t LASolverCholsky<T>::Inverse(const rsmat_t& S) const
{
    return inv(S);
}


template class LASolverEigen  <double>;
template class LASolverSVD    <double>;
template class LASolverCholsky<double>;
