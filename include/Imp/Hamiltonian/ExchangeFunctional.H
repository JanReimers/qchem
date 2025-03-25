// File: ExchangeFunctional.h  Interface for all ExchangeFunctional exchange functions.
#ifndef _ExchangeFunctional_H_
#define _ExchangeFunctional_H_



#include "Imp/Misc/pmstream.h"
#include <ScalarFunction.H>

class ExFunctional
    : public virtual PMStreamableObject
    , public virtual ScalarFunction<double>
{
public:
    virtual void            InsertChargeDensity(const Exact_CD*)      =0;
    virtual Vector<double>  GetVxcs(const Vector<double>& ChargeDensities) const=0;
    virtual double          GetVxc(                double ChargeDensity) const=0;
    virtual void            SetPolarized(bool)=0;

    virtual std::ostream&       Write(std::ostream&) const=0;
    virtual std::istream&       Read (std::istream&)      =0;
    static  ExFunctional* Factory(std::istream&);
};

#endif //_ExchangeFunctional_H_
