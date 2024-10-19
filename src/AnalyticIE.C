#include <AnalyticIE.H>
#include <LASolver.H>
#include <LAParams.H>
#include "oml/smatrix.h"

template <class T> typename AnalyticIE<T>::RSMat AnalyticIE<T>::
    MakeInverse(const RSMat& S,const LAParams& lap) 
{
    LASolver<double>* las=LASolver<double>::Factory(lap);
    RSMat Sinv=las->Inverse(S);
    delete las;
    return Sinv;
}

template class AnalyticIE<double>;
