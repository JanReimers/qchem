#include "IntegralEngine1.H"
#include "LASolver/LASolver.H"
#include "SCFIterator/IterationParams.H"

template <class T> typename IntegralEngine1<T>::RSMat IntegralEngine1<T>::
    MakeInverse(const RSMat& S) const
{
    LinearAlgebraParams lap={qchem::Lapack,qchem::SVD,1e-6,1e-12};
    LASolver<double>* las=LASolver<double>::Factory(lap);
    return las->Inverse(S);
//    return InvertSymmetric(S);
}

template class IntegralEngine1<double>;
