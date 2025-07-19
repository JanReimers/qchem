// File: IntegralConstraintFF.H  Integral constrained fit.
#ifndef _IntegralConstraintFF_H_
#define _IntegralConstraintFF_H_



#include "ConstrainedFF.H"

template <class T> class IntegralConstrainedFF
    : public ConstrainedFF<T>
{
public:
    typedef typename ConstrainedFF<T>::mesh_t mesh_t;
    typedef typename ConstrainedFF<T>::bs_t   bs_t;

    IntegralConstrainedFF(                                   );
    IntegralConstrainedFF(bs_t&, mesh_t& );
};

#endif //_IntegralConstraintFF_H_

