// File: EigenSolver.C  General eigen solver.

#include "Misc/EigenSolver.H"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include "oml/numeric.h"
#include "oml/cnumeric.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>

template <class T> EigenSolverOMLEigen<T>::EigenSolverOMLEigen(const SMat& S, double tolerance)
{
    V=Orthogonalize(S,tolerance);
}

template <class T> typename EigenSolver<T>::Mat 
EigenSolverOMLEigen<T>::Orthogonalize(const SMat& S, double tolerance) const
{
    auto [V,w] =Diagonalize(S);
//    RVec W(Wc.size());
//    for (unsigned int i=1;i<=Wc.size();i++) W(i)=real(Wc(i));
//    const RVec& w(W);
    cout << "Eigen values of S :" << w << std::endl;
//  cout << "Tolerance =" << Tolerance << std::endl;
//
//  Find out how many orht-normal linear combinations pass
//  The tolerance test.  Combinations with very small eigen values
//  Are linearly dependant.
//  I'm not sure this is true, I suspect that combinations with
//  very close eigen values are linearly dependant.
//
//  for (i=1;i<=n && W(i) < Tolerance; i++, NumAccepted--)
//    V.SetLimits(MatLimits(n,n));
//    V=S.SubMatrix(MatLimits(n,n));
////
//  Now rescale the columns in V V(i) = V(i)/sqrt(W(i)).
//
    size_t n=w.size();
    typename Mat::Subscriptor v(V);
    for (size_t i=1; i<=n; i++)
        for (size_t j=1; j<=n; j++)
        {
            if (w(i)<=0.0)
            {
                std::cerr << "Overlap matrix has a zero or negative eigen value " << w(i) << std::endl;
            }
            v(j,i)/=sqrt(w(i));
        }
//  cout << "V=" << V << std::endl;

    return V;
}

template <class T> typename EigenSolver<T>::UdType EigenSolverOMLEigen<T>::Solve(const SMat& Ham) const
{
    assert(!isnan(Ham));
	StreamableObject::SetToPretty();
	//std::cout << "Ham=" << Ham << std::endl;
	//std::cout << "HV=" << V << std::endl;
    Mat HPrime = ~V * Ham * V;  //Transform to orthogonal coordinates.
//	cout << "H'=" << HPrime << std::endl;
//	StreamableObject::SetToBinary();
    double del=MakeSymmetric(HPrime);
    SMat HS(HPrime.GetLimits());
    index_t rl=HS.GetLimits().Row.Low;
    for (index_t j:HS.cols())
        for (index_t i=rl;i<=j;i++)
            HS(i,j)=HPrime(i,j);
        

    if (fabs(del) > 1e-10)
    {
        std::cerr << "Warning: Hamiltonian asymmetry = " << del << " is big!" << std::endl;
//        exit (-1);
    }
//	cout << "HPrime=" << HPrime << std::endl;
    auto [U,e]  = Diagonalize(HS);        //Get eigen solution.
//    EigenValues.SetLimits(Wc.size());
//    for (unsigned int i=1;i<=Wc.size();i++) EigenValues(i)=real(Wc(i));
//	cout << "raw EigenVectors =" << HPrime << std::endl;
    U = V * U;                   //Back transform.
//	cout << "EigenVectors =" << EigenVectors << std::endl;
    return std::make_tuple(U,e);
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
//template class EigenSolver<std::complex<double> >;
