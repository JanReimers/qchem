// File: EigenSolverOML.C  General eigen solver.

#include "Misc/EigenSolver.H"
#include "Misc/EigenSolverOML.H"
#include "oml/smatrix.h"
#include "oml/numeric.h"
#include <iostream>

template <class T> typename EigenSolver<T>::UdType EigenSolverOMLCommon<T>::Solve(const SMat& Ham) const
{
    assert(!isnan(Ham));
	StreamableObject::SetToPretty();
    Mat HPrime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    double del=MakeSymmetric(HPrime);
    if (fabs(del) > 1e-10)
    {
        std::cerr << "Warning: Hamiltonian asymmetry = " << del << " is big!" << std::endl;
    }
    SMat HS(HPrime); //Make a symmetrix/Hermition version.
    auto [U,e]  = Diagonalize(HS);  //Get eigen solution.
    U = V * U;                      //Back transform.
    return std::make_tuple(U,e);
}


template <class T> EigenSolverOMLEigen<T>::EigenSolverOMLEigen(const SMat& S, double tolerance)
{
    auto [U,w] =Diagonalize(S);
    EigenSolverCommon<T>::Rescale(U,w);
    EigenSolverCommon<T>::AssignVs(U,~U);
}

template <class T> EigenSolverOMLSVD<T>::EigenSolverOMLSVD(const SMat& S, double tolerance)
{
    auto [U,s,V] =SVD(S);
//    Mat sM(s.GetLimits(),s.GetLimits());
//    Fill(sM,0.0);
//    sM.GetDiagonal()=s;
    //double err1=Max(fabs(U*sM*~V-S));
    EigenSolverCommon<T>::Rescale(U,s);
    EigenSolverCommon<T>::Rescale(V,s);
    //double err2=Max(fabs(U*~V-S));
    //std::cout << "SVD errors " << err1 << " " << err2 << std::endl;
    EigenSolverCommon<T>::AssignVs(U,~V);
}

template <class T> EigenSolverOMLCholsky<T>::EigenSolverOMLCholsky(const SMat& S, double tolerance)
{
    Mat U=S;
    Cholsky(U); //U is noe upper triangular, S=U*U_dagger 
    Mat Uinv=U; //Copy
    InvertTriangular(Uinv); //
//    double err1=Max(fabs(U*~U-S));
//    std::cout << "Cholsky errors " << err1 << std::endl;
    EigenSolverCommon<T>::AssignVs(~Uinv,Uinv);
}


template class EigenSolverOMLCommon<double>;
template class EigenSolverOMLEigen<double>;
template class EigenSolverOMLSVD<double>;
template class EigenSolverOMLCholsky<double>;
