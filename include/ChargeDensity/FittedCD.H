// File: FittedCD.H  Fitted charge density.
#ifndef _FittedCD_H_
#define _FittedCD_H_

#include <BasisSet/fwd.H>
#include <Fitting/FittedFunction.H>

//----------------------------------------------------------------------------------
//
//
class FittedCD
    : public virtual ScalarFunction<double>
    , public virtual FittedFunction
{
public:
    using FittedFunction::DoFit;
    typedef SMatrix<double> SMat;
    
    virtual double GetSelfRepulsion    (                       ) const=0;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual SMat   GetRepulsion(const TOrbital_DFT_IBS<double>*) const=0;
    //Required for creating a polarized CD from and un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
};


#endif //_FittedCD_H_
