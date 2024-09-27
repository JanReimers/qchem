// File: PolarizedCD.H  Interface for a spin polarized charge density.
#ifndef _PolarizedCD_H_
#define _PolarizedCD_H_

#include "ChargeDensity/ChargeDensity.H"

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//  All member functions just return the unpolarized answer.
//
class PolarizedCD
    : public virtual ChargeDensity
{
public:
    virtual       ChargeDensity* GetChargeDensity(const Spin&)      =0;
    virtual const ChargeDensity* GetChargeDensity(const Spin&) const=0;

    virtual SMat   GetOverlap  (const BasisSet*) const;
    virtual SMat   GetRepulsion(const BasisSet*) const;
    virtual SMat   GetExchange (const BasisSet*) const;

    virtual double GetEnergy(const HamiltonianTerm*) const;

    virtual double GetTotalCharge      (                     ) const;  // <ro>
    virtual double GetTotalSpin        (                     ) const;  // <up>-<down>

    virtual void   InjectOverlaps  (FittedFunction*, const BasisSet* theFitBasisSet) const;
    virtual void   InjectRepulsions(FittedFunction*, const BasisSet* theFitBasisSet) const;

    virtual bool   IsPolarized  (                           ) const;
    virtual void   ReScale      (double factor              )      ;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;  //Usefull for single atom charge densities.
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //Convergence check.

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;
    virtual void   Eval(const Mesh&, RVec&) const;
};

#endif //_PolarizedCD_H_
