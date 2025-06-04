// File: LASolverOML.C  General eigen solver.

#include <LASolver.H>
#include "Imp/LASolver/LASolverOML.H"
#include "oml/smatrix.h"
#include "oml/diagonalmatrix.h"
#include "oml/numeric.h"
#include <iostream>

template <class T> typename LASolver<T>::UdType LASolverOMLCommon<T>::Solve(const SMat& Ham) const
{
    assert(!isnan(Ham));
	Mat HPrime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    SMat HS=MakeSymmetric(HPrime,"Hamiltonian");
    auto [U,e]  = Diagonalize(HS);  //Get eigen solution.
    U = V * U;                      //Back transform.
    return std::make_tuple(U,e);
}

template <class T> typename LASolver<T>::UUdType LASolverOMLCommon<T>::SolveOrtho(const SMat& HPrime) const
{
    assert(!isnan(HPrime));
    auto [Uprime,e]  = Diagonalize(HPrime);  //Get eigen solution.
    Mat U = V * Uprime;                      //Back transform.
    return std::make_tuple(U,Uprime,e);
}


template <class T> void LASolverOMLEigen<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,w] =Diagonalize(S);
    LASolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    LASolverCommon<T>::Rescale(U,w);
    LASolverCommon<T>::AssignVs(U,~U);
    LASolverCommon<T>::Diag=w; //Preserve eigend values.
}

template <class T> typename LASolver<T>::RSMat LASolverOMLEigen<T>::Inverse(const RSMat& S) const
{
    auto [U,w] =Diagonalize(S);
    LASolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    DiagonalMatrix<T> winv(RVec(1.0/w));
    Mat Sfull(U*winv*~U);
    return MakeSymmetric(Sfull,"Inverse");
}


template <class T> void LASolverOMLSVD<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,s,V] =SVD(S);
    Mat Vt=~V;
    LASolverCommon<T>::Truncate(U,s,Vt,itsParams.TruncationTolerance);
    LASolverCommon<T>::Rescale(U,s,Vt);
    LASolverCommon<T>::AssignVs(U,Vt);
    LASolverCommon<T>::Diag=s; //Preserve SVs.
}

template <class T> typename LASolver<T>::RSMat LASolverOMLSVD<T>::Inverse(const RSMat& S) const
{
    auto [U,s,V] =SVD(S);
    Mat Vt=~V;
    LASolverCommon<T>::Truncate(U,s,Vt,itsParams.TruncationTolerance);
    DiagonalMatrix<T> sinv(RVec(1.0/s));
    Mat Sfull(~Vt*sinv*~U);
    return MakeSymmetric(Sfull,"Inverse");
}

template <class T> void LASolverOMLCholsky<T>::SetBasisOverlap(const SMat& S)
{
    Mat U=S;
    Cholsky(U); //U is noe upper triangular, S=U*U_dagger 
    Mat Uinv=U; //Copy
    InvertTriangular(Uinv); //
    LASolverCommon<T>::AssignVs(~Uinv,Uinv);
    LASolverCommon<T>::Diag=U.GetDiagonal();
}

template <class T> typename LASolver<T>::RSMat LASolverOMLCholsky<T>::Inverse(const RSMat& S) const
{
    return InvertSymmetric(S);
}



template class LASolverOMLCommon<double>;
template class LASolverOMLEigen<double>;
template class LASolverOMLSVD<double>;
template class LASolverOMLCholsky<double>;
