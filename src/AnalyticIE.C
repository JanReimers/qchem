#include "AnalyticIE.H"
#include "LASolver/LASolver.H"
#include "Imp/LAParams.H"
#include "oml/smatrix.h"

template <class T> typename AnalyticIE<T>::RSMat AnalyticIE<T>::
    MakeInverse(const RSMat& S,const LAParams& lap) 
{
    //LAParams lap={qchem::Lapack,qchem::SVD,1e-6,1e-12};
    LASolver<double>* las=LASolver<double>::Factory(lap);
    return las->Inverse(S);
}

template class AnalyticIE<double>;
