#include "AnalyticIE.H"
#include "LASolver/LASolver.H"
#include "SCFIterator/IterationParams.H"

template <class T> typename AnalyticIE<T>::RSMat AnalyticIE<T>::
    MakeInverse(const RSMat& S) const
{
    LinearAlgebraParams lap={qchem::Lapack,qchem::SVD,1e-6,1e-12};
    LASolver<double>* las=LASolver<double>::Factory(lap);
    return las->Inverse(S);
//    return InvertSymmetric(S);
}

template class AnalyticIE<double>;
