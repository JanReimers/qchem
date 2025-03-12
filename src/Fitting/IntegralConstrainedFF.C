// File: IntegralConstrainedFF.C  Integral constrained fit.



#include "Imp/Fitting/IntegralConstrainedFF.H"
#include <BasisSet.H>
#include "oml/matrix.h"

template <class T> IntegralConstrainedFF<T>::IntegralConstrainedFF()
    : ConstrainedFF<T>()
{};

template <class T> IntegralConstrainedFF<T>::
IntegralConstrainedFF(bs_t& theFitBasisSet, mesh_t&  m)
   // : ConstrainedFF<T>(theFitBasisSet,dynamic_cast<const TIrrepBasisSet<T>*>(theFitBasisSet.get())->GetDataBase()->GetCharge(theFitBasisSet.get()),m)
    : ConstrainedFF<T>(theFitBasisSet,theFitBasisSet->Charge(),m)
    {};

template class IntegralConstrainedFF<double>;
