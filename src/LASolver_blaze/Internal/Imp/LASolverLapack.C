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
template <class T> LASolverLapackCommon_blaze<T>::LASolverLapackCommon_blaze(const LAParams& lap)
    : LASolverCommon_blaze<T>(lap)
    {
    }
template <class T> LASolverLapackCommon_blaze<T>::~LASolverLapackCommon_blaze()
{
}

template <class T> typename LASolver_blaze<T>::Ud_t LASolverLapackCommon_blaze<T>::Solve(const smat_t& Ham) const
{
    assert(!isnan(Ham));
    mat_t Hprime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    smat_t Hsym=MakeSymmetric(Hprime,"Hamiltonian");
    rvec_t d;
    mat_t  U;
    blaze::eigen(Hsym,d,U);
    U = V * U;                      //Back transform.
    return std::make_tuple(U,d);
}
template <class T> typename LASolver_blaze<T>::UUd_t LASolverLapackCommon_blaze<T>::SolveOrtho(const smat_t& Hprime) const
{
    assert(!isnan(Hprime));
    rvec_t d;
    mat_t  Uprime;
    blaze::eigen(Hprime,d,Uprime);
    // auto [Uprime,e]  =itsLapackEigenSolver->SolveAll(HPrime,itsParams.abstol);  //Get eigen solution.
    mat_t U = V * Uprime;                      //Back transform.
    return std::make_tuple(U,Uprime,d);
}


template <class T> void LASolverLapackEigen_blaze<T>::SetBasisOverlap(const smat_t& S)
{
    rvec_t d;
    mat_t  U;
    blaze::eigen(S,d,U);
    // auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon_blaze<T>::Truncate(U,d,itsParams.TruncationTolerance);
    LASolverCommon_blaze<T>::Rescale(U,d);
    LASolverCommon_blaze<T>::AssignVs(U,trans(U));
    LASolverCommon_blaze<T>::Diag=d; //Preserve eigend values.
}

template <class T> typename LASolver_blaze<T>::rsmat_t LASolverLapackEigen_blaze<T>::Inverse(const rsmat_t& S) const
{
    rvec_t d;
    mat_t  U;
    blaze::eigen(S,d,U);
    // auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon_blaze<T>::Truncate(U,d,itsParams.TruncationTolerance);
    dmat_t winv(d.size());
    blaze::diagonal(winv)=1.0/d;
    mat_t Sfull(U*winv*trans(U));
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverLapackSVD_blaze<T>::SetBasisOverlap(const smat_t& S)
{
    rvec_t s;
    mat_t  U,Vt;
    blaze::svd(S,U,s,Vt);
    LASolverCommon_blaze<T>::Truncate(U,s,Vt,itsParams.TruncationTolerance);
    LASolverCommon_blaze<T>::Rescale(U,s,Vt);
    LASolverCommon_blaze<T>::AssignVs(U,Vt);
    LASolverCommon_blaze<T>::Diag=s; //Preserve SVs.
}

template <class T> typename LASolverLapackSVD_blaze<T>::rsmat_t LASolverLapackSVD_blaze<T>::Inverse(const rsmat_t& S) const
{
    rvec_t s;
    mat_t  U,Vt;
    blaze::svd(S,U,s,Vt);
    LASolverCommon_blaze<T>::Truncate(U,s,Vt,itsParams.TruncationTolerance);
    dmat_t sinv(s.size());
    blaze::diagonal(sinv)=1.0/s;
    mat_t  Sfull(trans(Vt)*sinv*trans(U));
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverLapackCholsky_blaze<T>::SetBasisOverlap(const smat_t& S)
{
    mat_t Sm(S);
    blaze::potrf( Sm, 'U' );
    blaze::trtri( Sm, 'U', 'N' );
    size_t N=Sm.rows();
    for (size_t i=0;i<N;i++)
        for (size_t j=i+1;j<N;j++)
            Sm(j,i)=0;
    umat_t V(Sm);
    

    LASolverCommon_blaze<T>::AssignVs(V,blaze::trans(V)); //V=U^-1, Vd=transpose(U^-1)
    LASolverCommon_blaze<T>::Diag=blaze::diagonal(V);
}

template <class T> typename LASolver_blaze<T>::rsmat_t LASolverLapackCholsky_blaze<T>::Inverse(const rsmat_t& S) const
{
    return inv(S);
}


template class LASolverLapackCommon_blaze <double>;
template class LASolverLapackEigen_blaze  <double>;
template class LASolverLapackSVD_blaze    <double>;
template class LASolverLapackCholsky_blaze<double>;
