// File: ExchangeFunctionalImplementation.H   Common implementation details for an exchange potential.
#ifndef _ExchangeFunctionalImplementation_H_
#define _ExchangeFunctionalImplementation_H_



#include "Hamiltonian/ExchangeFunctional.H"
#include "FunctionsImp/DumpScalarData.H"
#include "FunctionsImp/ScalarFunctionBuffer.H"

class ExchangeFunctionalImplementation
    : public  virtual ExchangeFunctional
    , private DumpScalarData
    , private ScalarFunctionBuffer<double>
{
public:
    ExchangeFunctionalImplementation(               );

    virtual void    InsertChargeDensity(const ChargeDensity*);
    virtual Vector<double> GetVxcs(const Vector<double>& ChargeDensities) const;
    virtual void            SetPolarized(bool p)
    {
        isPolarized=p;
    }

protected:
    virtual void  Eval    (const Mesh&, Vector<double>&) const;

    const ChargeDensity* itsChargeDensity;
    bool                 isPolarized;
};

#endif //_ExchangeFunctionalImplementation_H_
