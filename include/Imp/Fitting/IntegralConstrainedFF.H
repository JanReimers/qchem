// File: IntegralConstraintFF.H  Integral constrained fit.
#ifndef _IntegralConstraintFF_H_
#define _IntegralConstraintFF_H_



#include "Imp/Fitting/ConstrainedFF.H"

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    IntegralConstrainedFF(                                   );
    IntegralConstrainedFF(const rc_ptr<IrrepBasisSet>&, const rc_ptr<Mesh>& );
};

#endif //_IntegralConstraintFF_H_

