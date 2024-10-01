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
template <class T> EigenSolverLapackCommon<T>::EigenSolverLapackCommon()
    : itsLapackEigenSolver(new oml::LapackEigenSolver<T>())
    , itsLapackSVDSolver  (new oml::LapackSVDSolver  <T>())
    , eps(1e-12)
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
    auto [U,e]  =itsLapackEigenSolver->SolveAll(HPrime,eps);  //Get eigen solution.
    U = V * U;                      //Back transform.
    return std::make_tuple(U,e);
}


template <class T> EigenSolverLapackEigen<T>::EigenSolverLapackEigen(const SMat& S, double tolerance)
{
    auto [U,w] =itsLapackEigenSolver->SolveAll(Mat(S),eps);
    EigenSolverCommon<T>::Truncate(U,w,tolerance);
    EigenSolverCommon<T>::Rescale(U,w);
    EigenSolverCommon<T>::AssignVs(U,~U);
}

template <class T> EigenSolverLapackSVD<T>::EigenSolverLapackSVD(const SMat& S, double tolerance)
{
    auto [U,sM,Vt] =itsLapackSVDSolver->SolveAll(S,eps);
    EigenSolverCommon<T>::Truncate(U,sM,Vt,tolerance);
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

template <class T> EigenSolverLapackCholsky<T>::EigenSolverLapackCholsky(const SMat& S, double tolerance)
{
    std::cerr << "General Eigen solver Lapack Cholsky is not implemented yet" << std::endl;
    exit(-1);
}


template class EigenSolverLapackEigen  <double>;
template class EigenSolverLapackSVD    <double>;
template class EigenSolverLapackCholsky<double>;
