// File: LASolver_blaze/Internal/Imp/LASolverLapack.C  General eigen solver.
module;
#include <cassert>
#include <iostream>
#include "blaze/Math.h" 
module qchem.LASolver_blaze.Internal.Lapack;

//-----------------------------------------------------------------------------
//
//  Lapack specific code
//
template <class T> LASolverLapackCommon_blaze<T>::LASolverLapackCommon_blaze(double truncationTolerance)
    : LASolverCommon<T>(truncationTolerance)
    {
    }
template <class T> LASolverLapackCommon_blaze<T>::~LASolverLapackCommon_blaze()
{
}

template <class T> typename LASolver<T>::Ud_t LASolverLapackCommon_blaze<T>::Solve(const smat_t<T>& Ham) const
{
    assert(!isnan(Ham));
    mat_t<T> Hprime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    smat_t<T> Hsym=MakeSymmetric(Hprime,"Hamiltonian");
    rvec_t d;
    mat_t<T>  U;
    blaze::eigen(Hsym,d,U);
    U = V * U;                      //Back transform.
    return std::make_tuple(U,d);
}
template <class T> typename LASolver<T>::UUd_t LASolverLapackCommon_blaze<T>::SolveOrtho(const smat_t<T>& Hprime) const
{
    assert(!isnan(Hprime));
    rvec_t d;
    mat_t<T>  Uprime;
    blaze::eigen(Hprime,d,Uprime);
    // auto [Uprime,e]  =itsLapackEigenSolver->SolveAll(HPrime,itsParams.abstol);  //Get eigen solution.
    mat_t<T> U = V * Uprime;                      //Back transform.
    return std::make_tuple(U,Uprime,d);
}


template <class T> void LASolverLapackEigen_blaze<T>::SetBasisOverlap(const smat_t<T>& S)
{
    rvec_t d;
    mat_t<T>  U;
    blaze::eigen(S,d,U);
    // auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon<T>::Truncate(U,d,itsTruncationTolerance);
    LASolverCommon<T>::Rescale(U,d);
    LASolverCommon<T>::AssignVs(U,trans(U));
    LASolverCommon<T>::Diag=d; //Preserve eigend values.
}

template <class T> rsmat_t LASolverLapackEigen_blaze<T>::Inverse(const rsmat_t& S) const
{
    rvec_t d;
    mat_t<T>  U;
    blaze::eigen(S,d,U);
    // auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon<T>::Truncate(U,d,itsTruncationTolerance);
    blaze::DiagonalMatrix<rmat_t> winv(d.size());
    blaze::diagonal(winv)=1.0/d;
    mat_t<T> Sfull(U*winv*trans(U));
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverLapackSVD_blaze<T>::SetBasisOverlap(const smat_t<T>& S)
{
    rvec_t s;
    mat_t<T>  U,Vt;
    blaze::svd(S,U,s,Vt);
    LASolverCommon<T>::Truncate(U,s,Vt,itsTruncationTolerance);
    LASolverCommon<T>::Rescale(U,s,Vt);
    LASolverCommon<T>::AssignVs(U,Vt);
    LASolverCommon<T>::Diag=s; //Preserve SVs.
}

template <class T> rsmat_t LASolverLapackSVD_blaze<T>::Inverse(const rsmat_t& S) const
{
    rvec_t s;
    mat_t<T>  U,Vt;
    blaze::svd(S,U,s,Vt);
    LASolverCommon<T>::Truncate(U,s,Vt,itsTruncationTolerance);
    blaze::DiagonalMatrix<rmat_t> sinv(s.size());
    blaze::diagonal(sinv)=1.0/s;
    mat_t<T>  Sfull(trans(Vt)*sinv*trans(U));
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverLapackCholsky_blaze<T>::SetBasisOverlap(const smat_t<T>& S)
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

template <class T> rsmat_t LASolverLapackCholsky_blaze<T>::Inverse(const rsmat_t& S) const
{
    return inv(S);
}


template class LASolverLapackCommon_blaze <double>;
template class LASolverLapackEigen_blaze  <double>;
template class LASolverLapackSVD_blaze    <double>;
template class LASolverLapackCholsky_blaze<double>;
