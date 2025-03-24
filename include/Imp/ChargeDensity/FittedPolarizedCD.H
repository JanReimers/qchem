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
    FittedPolarizedCD();
    FittedPolarizedCD(const FittedCD* unpolcd, double Stotal);
    FittedPolarizedCD(FittedCD* up, FittedCD* down);
    FittedPolarizedCD(const FittedPolarizedCD&);
    ~FittedPolarizedCD();

    virtual       FittedCD* GetChargeDensity(const Spin&)      ;
    virtual const FittedCD* GetChargeDensity(const Spin&) const;

    virtual double GetSelfRepulsion    (                     ) const;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual double GetRepulsion        (const FittedFunction*) const;  // 1/2 <ro(1) | 1/r12 | ff(2)>
    virtual Vector<double> GetRepulsion3C(const Fit_IBS* theFitBasisSet) const;
    virtual SMat   GetRepulsion(const TOrbital_DFT_IBS<double>*) const;

    virtual double DoFit           (const ScalarFFClient&);
    virtual double DoFit           (const DensityFFClient&);
    virtual void   ReScale         (double factor      );  //Ro *= factor
    virtual void   ShiftOrigin     (const RVec3&       );  //Usefull for single atom charge densities.
    virtual void   FitMixIn        (const FittedFunction&,double)      ; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const;

    virtual FittedCD* Clone(        ) const;

private:
    FittedCD* itsSpinUpCD;
    FittedCD* itsSpinDownCD;
};



#endif // _FittedPolarizedCD_H_
