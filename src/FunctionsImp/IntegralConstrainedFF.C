// File: IntegralConstrainedFF.C  Integral constrained fit.



#include "FunctionsImp/IntegralConstrainedFF.H"
#include "BasisSet.H"
#include "IntegralDataBase.H"
#include "oml/matrix.h"

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(const rc_ptr<IrrepBasisSet>& theFitBasisSet, bool CDfit)
    : ConstrainedFF<T>(theFitBasisSet,dynamic_cast<const TIrrepBasisSet<T>*>(theFitBasisSet.get())->GetDataBase()->GetCharge(theFitBasisSet.get()),CDfit)
{};

template class IntegralConstrainedFF<double>;
