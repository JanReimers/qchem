// File: ExchangeFunctionalImplementation.H   Common implementation details for an exchange potential.
#ifndef _ExFunctionalImp_H_
#define _ExFunctionalImp_H_



#include "ExchangeFunctional.H"

class ExFunctionalImp
    : public  virtual ExFunctional
{
public:
    ExFunctionalImp(               );

    virtual void           InsertChargeDensity(const DM_CD*);
    virtual Vector<double> GetVxcs(const Vector<double>& ChargeDensities) const;
    virtual void           SetPolarized(bool p) {isPolarized=p;}

protected:
    // virtual void  Eval    (const Mesh&, Vector<double>&) const;

    const DM_CD* itsChargeDensity;
    bool            isPolarized;
};

#endif //_ExFunctionalImp_H_
