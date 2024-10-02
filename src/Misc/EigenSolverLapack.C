// File: EigenSolverLapack.C  General eigen solver.

#include "Misc/EigenSolver.H"
#include "Misc/EigenSolverLapack.H"
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
template <class T> EigenSolverLapackCommon<T>::EigenSolverLapackCommon(const LinearAlgebraParams& lap)
    : EigenSolverCommon<T>(lap)
    , itsLapackEigenSolver(new oml::LapackEigenSolver<T>())
    , itsLapackSVDSolver  (new oml::LapackSVDSolver  <T>())
    {
        assert(itsLapackEigenSolver);
        assert(itsLapackSVDSolver);
    }


template <class T> typename EigenSolver<T>::UdType EigenSolverLapackCommon<T>::Solve(const SMat& Ham) const
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


template <class T> void EigenSolverLapackEigen<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),itsParams.abstol);
    EigenSolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    EigenSolverCommon<T>::Rescale(U,w);
    EigenSolverCommon<T>::AssignVs(U,~U);
}

template <class T> void EigenSolverLapackSVD<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,sM,Vt] =itsLapackSVDSolver->SolveAll(S,itsParams.abstol);
    EigenSolverCommon<T>::Truncate(U,sM,Vt,itsParams.TruncationTolerance);
    RVec s=sM.GetDiagonal();
//    Mat sM(s.GetLimits(),s.GetLimits());
//    Fill(sM,0.0);
//    sM.GetDiagonal()=s;
    //double err1=Max(fabs(U*sM*~V-S));
    EigenSolverCommon<T>::Rescale(U,s,Vt);
    //double err2=Max(fabs(U*~V-S));
    //std::cout << "SVD errors " << err1 << " " << err2 << std::endl;
    EigenSolverCommon<T>::AssignVs(U,Vt);
}

template <class T> void EigenSolverLapackCholsky<T>::SetBasisOverlap(const SMat& S)
{
    std::cerr << "General Eigen solver Lapack Cholsky is not implemented yet" << std::endl;
    exit(-1);
}


template class EigenSolverLapackEigen  <double>;
template class EigenSolverLapackSVD    <double>;
template class EigenSolverLapackCholsky<double>;
