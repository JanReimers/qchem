// File: LASolverOML.C  General eigen solver.

#include "LASolver/LASolver.H"
#include "Imp/LASolver/LASolverOML.H"
#include "oml/smatrix.h"
#include "oml/numeric.h"
#include <iostream>

template <class T> typename LASolver<T>::UdType LASolverOMLCommon<T>::Solve(const SMat& Ham) const
{
    assert(!isnan(Ham));
	StreamableObject::SetToPretty();
    Mat HPrime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    double del=MakeSymmetric(HPrime);
//    if (fabs(del) > 1e-10)
//    {
//        std::cerr << "Warning: Hamiltonian asymmetry = " << del << " is big!" << std::endl;
//    }
    SMat HS(HPrime); //Make a symmetrix/Hermition version.
    auto [U,e]  = Diagonalize(HS);  //Get eigen solution.
    U = V * U;                      //Back transform.
    return std::make_tuple(U,e);
}


template <class T> void LASolverOMLEigen<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,w] =Diagonalize(S);
    LASolverCommon<T>::Truncate(U,w,itsParams.TruncationTolerance);
    LASolverCommon<T>::Rescale(U,w);
    LASolverCommon<T>::AssignVs(U,~U);
}

template <class T> void LASolverOMLSVD<T>::SetBasisOverlap(const SMat& S)
{
    auto [U,s,V] =SVD(S);
    LASolverCommon<T>::Truncate(U,s,V,itsParams.TruncationTolerance);

//    Mat sM(s.GetLimits(),s.GetLimits());
//    Fill(sM,0.0);
//    sM.GetDiagonal()=s;
    //double err1=Max(fabs(U*sM*~V-S));
    LASolverCommon<T>::Rescale(U,s);
    LASolverCommon<T>::Rescale(V,s);
    //double err2=Max(fabs(U*~V-S));
    //std::cout << "SVD errors " << err1 << " " << err2 << std::endl;
    LASolverCommon<T>::AssignVs(U,~V);
}

template <class T> void LASolverOMLCholsky<T>::SetBasisOverlap(const SMat& S)
{
    Mat U=S;
    Cholsky(U); //U is noe upper triangular, S=U*U_dagger 
    Mat Uinv=U; //Copy
    InvertTriangular(Uinv); //
//    double err1=Max(fabs(U*~U-S));
//    std::cout << "Cholsky errors " << err1 << std::endl;
    LASolverCommon<T>::AssignVs(~Uinv,Uinv);
}


template class LASolverOMLCommon<double>;
template class LASolverOMLEigen<double>;
template class LASolverOMLSVD<double>;
template class LASolverOMLCholsky<double>;
