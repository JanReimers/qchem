// File: IntegralConstrainedFF.C  Integral constrained fit.



#include "IntegralConstrainedFF.H"
#include <BasisSet/Fit_IBS.H>
#include "oml/matrix.h"

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(bs_t& fbs, mesh_t&  m)
    : ConstrainedFF<T>(fbs,fbs->Charge(),m)
    {};

template class IntegralConstrainedFF<double>;
