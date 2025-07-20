// File: LASolverLapack.C  General eigen solver.
module;
#include <cassert>
#include <iostream>
module qchem.LASolver.Internal.Lapack;

//-----------------------------------------------------------------------------
//
//  Lapack specific code
//
template <class T> LASolverLapackCommon<T>::LASolverLapackCommon(const LAParams& lap)
    : LASolverCommon<T>(lap)
    , itsLapackEigenSolver(new oml::LapackEigenSolver<T>())
    , itsLapackSVDSolver  (new oml::LapackSVDSolver  <T>())
    {
        assert(itsLapackEigenSolver);
        assert(itsLapackSVDSolver);
    }
template <class T> LASolverLapackCommon<T>::~LASolverLapackCommon()
{
    delete itsLapackEigenSolver;
    delete itsLapackSVDSolver;
}

template <class T> typename LASolver<T>::UdType LASolverLapackCommon<T>::Solve(const SMat& Ham) const
{
    assert(!isnan(Ham));
    Mat HPrime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    MakeSymmetric(HPrime,"Hamiltonian");
    auto [U,e]  =itsLapackEigenSolver->SolveAll(HPrime,itsParams.abstol);  //Get eigen solution.
    U = V * U;                      //Back transform.
    return std::make_tuple(U,e);
}
template <class T> typename LASolver<T>::UUdType LASolverLapackCommon<T>::SolveOrtho(const SMat& HPrime) const
{
    assert(!isnan(HPrime));
    auto [Uprime,e]  =itsLapackEigenSolver->SolveAll(HPrime,itsParams.abstol);  //Get eigen solution.
    Mat U = V * Uprime;                      //Back transform.
    return std::make_tuple(U,Uprime,e);
}


template <class T> void LASolverLapackEigen<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    LASolverCommon<T>::Rescale(U,w);
    LASolverCommon<T>::AssignVs(U,~U);
    LASolverCommon<T>::Diag=w; //Preserve eigend values.
}

template <class T> typename LASolver<T>::RSMat LASolverLapackEigen<T>::Inverse(const RSMat& S) const
{
    auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    DiagonalMatrix<T> winv(RVec(1.0/w));
    Mat Sfull(U*winv*~U);
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverLapackSVD<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,sM,Vt] =itsLapackSVDSolver->SolveAll(S);
    RVec s=sM.GetDiagonal();
    LASolverCommon<T>::Truncate(U,s,Vt,itsParams.TruncationTolerance);
    LASolverCommon<T>::Rescale(U,s,Vt);
    LASolverCommon<T>::AssignVs(U,Vt);
    LASolverCommon<T>::Diag=s; //Preserve SVs.
}

template <class T> typename LASolverLapackSVD<T>::RSMat LASolverLapackSVD<T>::Inverse(const RSMat& S) const
{
    auto [U,sM,Vt] =itsLapackSVDSolver->SolveAll(S);
    RVec s=sM.GetDiagonal();
    LASolverCommon<T>::Truncate(U,s,Vt,itsParams.TruncationTolerance);
    DiagonalMatrix<T> sinv(RVec(1.0/s));
    Mat Sfull(~Vt*sinv*~U);
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverLapackCholsky<T>::SetBasisOverlap(const SMat& S)
{
    Mat U=oml::LapackCholsky(S); // S=Ud*U
    Mat Uinv=oml::LapackInvertTriangular(U); //V=U^-1
    LASolverCommon<T>::AssignVs(Uinv,~Uinv); //V=U^-1, Vd=transpose(U^-1)
    LASolverCommon<T>::Diag=U.GetDiagonal();
}

template <class T> typename LASolver<T>::RSMat LASolverLapackCholsky<T>::Inverse(const RSMat& S) const
{
    return oml::LapackInvertSymmetric(S);
}


template class LASolverLapackCommon <double>;
template class LASolverLapackEigen  <double>;
template class LASolverLapackSVD    <double>;
template class LASolverLapackCholsky<double>;
