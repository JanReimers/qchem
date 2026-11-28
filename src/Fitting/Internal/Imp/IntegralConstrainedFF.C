// File: IntegralConstrainedFF.C  Integral constrained fit.
module;
module qchem.Fitting.Internal.FunctionFitterImp;
import qchem.Fitting.Types;

namespace qchem::Fitting
{

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(fbs_t& fbs, bool constrained)
    : ConstrainedFF<T>(fbs,fbs->Charge(),constrained)
    {};

template class IntegralConstrainedFF<double>;

} //namespace
