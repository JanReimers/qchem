// File: FittedPolarizedCD.H Implementation
#ifndef _FittedPolarizedCD_H_
#define _FittedPolarizedCD_H_

#include <ChargeDensity.H>

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//
class FittedPolarizedCD
    : public virtual FittedCD
    , public virtual Polarized_CD
    
{
public:
    FittedPolarizedCD(); // No UT coverage
    FittedPolarizedCD(const FittedCD* unpolcd, double Stotal); // No UT coverage
    FittedPolarizedCD(FittedCD* up, FittedCD* down); // No UT coverage
    FittedPolarizedCD(const FittedPolarizedCD&); // No UT coverage
    ~FittedPolarizedCD();

    virtual       FittedCD* GetChargeDensity(const Spin&)      ; // No UT coverage
    virtual const FittedCD* GetChargeDensity(const Spin&) const; // No UT coverage

    virtual double GetSelfRepulsion    (                     ) const;  // No UT coverage // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual double GetRepulsion        (const FittedFunction*) const;  // No UT coverage// 1/2 <ro(1) | 1/r12 | ff(2)>
    virtual Vector<double> GetRepulsion3C(const Fit_IBS* theFitBasisSet) const; // No UT coverage
    virtual SMat   GetRepulsion(const TOrbital_DFT_IBS<double>*) const; // No UT coverage

    virtual double DoFit           (const ScalarFFClient&); // No UT coverage
    virtual double DoFit           (const DensityFFClient&); // No UT coverage
    virtual void   ReScale         (double factor      );  // No UT coverage//Ro *= factor
    virtual void   ShiftOrigin     (const RVec3&       );  // No UT coverage//Usefull for single atom charge densities.
    virtual void   FitMixIn        (const FittedFunction&,double)      ; // No UT coverage// this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const; // No UT coverage

    virtual FittedCD* Clone(        ) const;

private:
    FittedCD* itsSpinUpCD;
    FittedCD* itsSpinDownCD;
};



#endif // _FittedPolarizedCD_H_
