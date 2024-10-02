// File: LASolverLapack.C  General eigen solver.

#include "LASolver/LASolver.H"
#include "Imp/LASolver/LASolverLapack.H"
#include "oml/numeric/LapackEigenSolver.H"
#include "oml/numeric/LapackSVDSolver.H"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include "oml/diagonalmatrix.h"
#include "oml/vector.h"
//#include <iostream>


//-----------------------------------------------------------------------------
//
//  Lapack specific code
//
template <class T> LASolverLapackCommon<T>::LASolverLapackCommon(const LinearAlgebraParams& lap)
    : LASolverCommon<T>(lap)
    , itsLapackEigenSolver(new oml::LapackEigenSolver<T>())
    , itsLapackSVDSolver  (new oml::LapackSVDSolver  <T>())
    {
        assert(itsLapackLASolver);
        assert(itsLapackSVDSolver);
    }


template <class T> typename LASolver<T>::UdType LASolverLapackCommon<T>::Solve(const SMat& Ham) const
{
    assert(!isnan(Ham));
	StreamableObject::SetToPretty();
    Mat HPrime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    double del=MakeSymmetric(HPrime);
//    if (fabs(del) > 1e-10)
//    {
//        std::cerr << "Warning: Hamiltonian asymmetry = " << del << " is big!" << std::endl;
//    }
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

template <class T> void LASolverLapackSVD<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,sM,Vt] =itsLapackSVDSolver->SolveAll(S,itsParams.abstol);
    LASolverCommon<T>::Truncate(U,sM,Vt,itsParams.TruncationTolerance);
    RVec s=sM.GetDiagonal();
//    Mat sM(s.GetLimits(),s.GetLimits());
//    Fill(sM,0.0);
//    sM.GetDiagonal()=s;
    //double err1=Max(fabs(U*sM*~V-S));
    LASolverCommon<T>::Rescale(U,s,Vt);
    //double err2=Max(fabs(U*~V-S));
    //std::cout << "SVD errors " << err1 << " " << err2 << std::endl;
    LASolverCommon<T>::AssignVs(U,Vt);
}

template <class T> void LASolverLapackCholsky<T>::SetBasisOverlap(const SMat& S)
{
    std::cerr << "General Eigen solver Lapack Cholsky is not implemented yet" << std::endl;
    exit(-1);
}


template class LASolverLapackEigen  <double>;
template class LASolverLapackSVD    <double>;
template class LASolverLapackCholsky<double>;
