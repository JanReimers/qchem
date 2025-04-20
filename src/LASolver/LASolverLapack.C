// File: LASolverLapack.C  General eigen solver.

#include <LASolver.H>
#include "Imp/LASolver/LASolverLapack.H"
#include "oml/numeric/LapackEigenSolver.H"
#include "oml/numeric/LapackSVDSolver.H"
#include "oml/numeric/LapackCholsky.H"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include "oml/diagonalmatrix.h"
#include "oml/vector.h"
//#include <iostream>


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
	StreamableObject::SetToPretty();
    Mat HPrime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    MakeSymmetric(HPrime,"Hamiltonian");
    auto [U,e]  =itsLapackEigenSolver->SolveAll(HPrime,itsParams.abstol);  //Get eigen solution.
    U = V * U;                      //Back transform.
    return std::make_tuple(U,e);
}


template <class T> void LASolverLapackEigen<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    LASolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    LASolverCommon<T>::Rescale(U,w);
    LASolverCommon<T>::AssignVs(U,~U);
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
    Mat U=oml::LapackCholsky(S);
    Mat Uinv=oml::LapackInvertTriangular(U); 
    LASolverCommon<T>::AssignVs(Uinv,~Uinv);
}

template <class T> typename LASolver<T>::RSMat LASolverLapackCholsky<T>::Inverse(const RSMat& S) const
{
    return oml::LapackInvertSymmetric(S);
}


template class LASolverLapackEigen  <double>;
template class LASolverLapackSVD    <double>;
template class LASolverLapackCholsky<double>;
