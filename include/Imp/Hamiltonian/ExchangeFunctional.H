// File: ExchangeFunctional.h  Interface for all ExchangeFunctional exchange functions.
#ifndef _ExchangeFunctional_H_
#define _ExchangeFunctional_H_


#include <iosfwd>
import qchem.ChargeDensity;
import qchem.ScalarFunction;
import qchem.Streamable;

class ExFunctional
    : public virtual Streamable
    , public virtual ScalarFunction<double>
{
public:
    virtual void            InsertChargeDensity(const DM_CD*)      =0;
    virtual Vector<double>  GetVxcs(const Vector<double>& ChargeDensities) const=0;
    virtual double          GetVxc(                double ChargeDensity) const=0;
    virtual void            SetPolarized(bool)=0;

    virtual std::ostream&       Write(std::ostream&) const=0;
};

#endif //_ExchangeFunctional_H_
