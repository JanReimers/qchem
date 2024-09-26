// File: EigenSolver.C  General eigen solver.



#include "Misc/EigenSolver.H"
#include "BasisSet/TBasisSet.H"
#include "BasisSet/IntegralDataBase.H"
//#include "oml/isnan.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include "oml/numeric.h"
#include "oml/cnumeric.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>

//inline double real(double f) {return f;}

template <class T> EigenSolver<T>::EigenSolver(const TBasisSet<T>& set, double Tolerance)
{
    Tolerance=0;
    index_t n=set.GetNumFunctions(), i;
//
//  Make basis set overlap matrix.
//
    IntegralDataBase<T>* db=set.GetDataBase();
    Mat S=db->GetOverlap();
//  StreamableObject::SetToPretty();
//  cout << "Normalized basis set overlap matrix :" << S;
//
//  Get eigen values and eigen vectors (now in S).
//
    Vec Wc=Diagonalize(S);
    RVec W(Wc.size());
    for (int i=1;i<=Wc.size();i++) W(i)=real(Wc(i));
    const RVec& w(W);
//  cout << "Eigen values of S :" << W << std::endl;
//  cout << "Tolerance =" << Tolerance << std::endl;
//
//  Find out how many orht-normal linear combinations pass
//  The tolerance test.  Combinations with very small eigen values
//  Are linearly dependant.
//  I'm not sure this is true, I suspect that combinations with
//  very close eigen values are linearly dependant.
//
//  for (i=1;i<=n && W(i) < Tolerance; i++, NumAccepted--)
    V.SetLimits(MatLimits(n,n));
    V=S.SubMatrix(MatLimits(n,n));
//
//  Now rescale the columns in V V(i) = V(i)/sqrt(W(i)).
//
    typename Mat::Subscriptor v(V);
    for (i=1; i<=n; i++)
        for (index_t j=1; j<=n; j++)
        {
            if (w(i)<=0.0)
            {
                std::cerr << "Overlap matrix has a zero or negative eigen value " << w(i) << std::endl;
            }
            v(j,i)/=sqrt(w(i));
        }
//  cout << "V=" << V << std::endl;

//
//  Now set the sizes of the eigen vectors and eigenvalues
//
    EigenVectors.SetLimits(MatLimits(n,n));
    EigenValues .SetLimits(VecLimits(  n));
}

template <class T> void EigenSolver<T>::Solve(const SMat& Ham)
{
    assert(!isnan(Ham));
	StreamableObject::SetToPretty();
	//std::cout << "Ham=" << Ham << std::endl;
	//std::cout << "HV=" << V << std::endl;
    Mat HPrime = ~V * Ham * V;  //Transform to orthogonal coordinates.
//	cout << "H'=" << HPrime << std::endl;
//	StreamableObject::SetToBinary();
    double del=MakeSymmetric(HPrime);

    if (fabs(del) > 1e-10)
    {
        std::cerr << "Warning: Hamiltonian asymmetry = " << del << " is big!" << std::endl;
//        exit (-1);
    }
//	cout << "HPrime=" << HPrime << std::endl;
    Vec Wc  = Diagonalize(HPrime);        //Get eigen solution.
    EigenValues.SetLimits(Wc.size());
    for (int i=1;i<=Wc.size();i++) EigenValues(i)=real(Wc(i));
//	cout << "raw EigenVectors =" << HPrime << std::endl;
    EigenVectors = V * HPrime;                   //Back transform.
//	cout << "EigenVectors =" << EigenVectors << std::endl;
}

template <class T> const typename EigenSolver<T>::Mat& EigenSolver<T>::GetEigenVectors() const
{
    return EigenVectors;
}

template <class T>  const typename EigenSolver<T>::RVec& EigenSolver<T>::GetEigenValues() const
{
    return EigenValues;
}

template class EigenSolver<double>;
//template class EigenSolver<std::complex<double> >;
