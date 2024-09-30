// File: EigenSolver.C  General eigen solver.

#include "Misc/EigenSolver.H"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include "oml/numeric.h"
#include "oml/cnumeric.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>

//-------------------------------------------------------------------------
//
//  Common level.
//
template <class T> void EigenSolverCommon<T>::Rescale(Mat& V,const RVec& w)
{
    for (auto j:V.cols())
        V.GetColumn(j)/=sqrt(w(j));
        
}

//-----------------------------------------------------------------------------
//
//  OML specific code
//
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
    Mat sM(s.GetLimits(),s.GetLimits());
    Fill(sM,0.0);
    sM.GetDiagonal()=s;
    //double err1=Max(fabs(U*sM*~V-S));
    EigenSolverCommon<T>::Rescale(U,s);
    EigenSolverCommon<T>::Rescale(V,s);
    //double err2=Max(fabs(U*~V-S));
    //std::cout << "SVD errors " << err1 << " " << err2 << std::endl;
    EigenSolverCommon<T>::AssignVs(U,~V);
}


template <class T> EigenSolver<T>* EigenSolver<T>::
    Factory(EigenSolver<T>::Pkg pkg,EigenSolver<T>::Ortho ortho,const SMat& S, double tolerance)
{
    EigenSolver<T>* ret=0;
    switch (pkg)
    {
    case OML:
        switch (ortho)
        {
        case Cholsky :
            break;
        case Eigen :
            ret=new EigenSolverOMLEigen<T>(S,tolerance);
            break;
        case SVD :
            ret=new EigenSolverOMLSVD<T>(S,tolerance);
            break;
        }
        break;
   case Lapack:
        switch (ortho)
        {
        case Cholsky :
            break;
        case Eigen :
            break;
        case SVD :
            break;
        break;
        }
    }
    return ret;
}

template class EigenSolver<double>;
template class EigenSolverCommon<double>;
template class EigenSolverOMLEigen<double>;
template class EigenSolverOMLSVD<double>;
//template class EigenSolver<std::complex<double> >;
