// File: BasisFunction.H  Abstract basis function interface.
#ifndef _BasisFunction_H_
#define _BasisFunction_H_



#include "Misc/Persistent/Persistent.H"

//--------------------------------------------------------------
//
//  This interface of the basis function provides no notion
//  of weather the function is real or complex (or vector!)
//  valued.  Normalization and total charge should always
//  be real.
//
class BasisFunction
    : public Persistent
{
public:
    virtual ~BasisFunction()  {};

    virtual bool   operator==(const BasisFunction&) const=0;

    virtual double GetNormalization() const=0;
    virtual double GetCharge       () const=0;

    virtual BasisFunction* Clone  (             ) const=0;
    static  BasisFunction* Factory(std::istream&)        ;
};

#endif //_BasisFunction_H_
