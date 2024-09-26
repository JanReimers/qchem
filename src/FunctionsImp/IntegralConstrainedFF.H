// File: IntegralConstraintFF.H  Integral constrained fit.
#ifndef _IntegralConstraintFF_H_
#define _IntegralConstraintFF_H_



#include "FunctionsImp/ConstrainedFF.H"

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    IntegralConstrainedFF(                                        );
    IntegralConstrainedFF(const rc_ptr<BasisSet>&,bool CDfit=false);
};

#endif //_IntegralConstraintFF_H_

