// File: IntegralConstrainedFF.C  Integral constrained fit.



#include "FunctionsImp/IntegralConstrainedFF.H"
#include "BasisSet/BasisSet.H"
#include "BasisSet/IntegralDataBase.H"
#include "oml/matrix.h"

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(const rc_ptr<BasisSet>& theFitBasisSet, bool CDfit)
    : ConstrainedFF<T>(theFitBasisSet,dynamic_cast<const TBasisSet<T>*>(theFitBasisSet.get())->GetDataBase()->GetCharge(),CDfit)
{};

template class IntegralConstrainedFF<double>;
