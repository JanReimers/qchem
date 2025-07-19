// File: IntegralConstrainedFF.C  Integral constrained fit.


#include <memory>
#include <vector>

#include "IntegralConstrainedFF.H"
import qchem.Fit_IBS;

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(bs_t& fbs, mesh_t&  m)
    : ConstrainedFF<T>(fbs,fbs->Charge(),m)
    {};

template class IntegralConstrainedFF<double>;
