// File: IntegralConstrainedFF.C  Integral constrained fit.
module;
module qchem.Fitting.Internal.FunctionFitterImp;
import qchem.Fitting.Types;

namespace qchem::Fitting
{

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(bs_t& fbs)
    : ConstrainedFF<T>(fbs,fbs->Charge())
    {};

template class IntegralConstrainedFF<double>;

} //namespace