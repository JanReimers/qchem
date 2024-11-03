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
    , public virtual PolarizedCD
    
{
public:
    FittedPolarizedCD();
    FittedPolarizedCD(const ChargeDensity* unpolcd, double Stotal);
    FittedPolarizedCD(ChargeDensity* up, ChargeDensity* down);
    FittedPolarizedCD(const FittedPolarizedCD&);
    ~FittedPolarizedCD();

    virtual       ChargeDensity* GetChargeDensity(const Spin&)      ;
    virtual const ChargeDensity* GetChargeDensity(const Spin&) const;

    virtual double GetSelfRepulsion    (                     ) const;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual double GetRepulsion        (const FittedFunction*) const;  // 1/2 <ro(1) | 1/r12 | ff(2)>

    virtual double DoFit           (const ScalarFFClient&);
    virtual double DoFit           (const DensityFFClient&);
    virtual void   ReScale         (double factor      );  //Ro *= factor
    virtual void   ShiftOrigin     (const RVec3&       );  //Usefull for single atom charge densities.
    virtual void   FitMixIn        (const FittedFunction&,double)      ; // this = this*(1-c) + that*c.
    virtual double FitGetChangeFrom(const FittedFunction&       ) const;

    virtual Vector<double> GetRepulsion3C(const IrrepBasisSet* theFitBasisSet) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual FittedCD* Clone(        ) const;

private:
    FittedCD* itsSpinUpCD;
    FittedCD* itsSpinDownCD;
};



#endif // _FittedPolarizedCD_H_
