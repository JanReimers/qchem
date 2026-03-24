// File: IntegralConstrainedFF.C  Integral constrained fit.
module;
#include <memory>
#include <vector>
module qchem.FittedFunctionImp;
import qchem.Fit_IBS;
import qchem.Conversions;

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(bs_t& fbs, mesh_t&  m)
    : ConstrainedFF<T>(fbs,convert(fbs->Charge()),m)
    {};

template class IntegralConstrainedFF<double>;
