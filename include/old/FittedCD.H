// File: FittedCD.H  Interface for the charge density category.
#ifndef _FittedCD_H_
#define _FittedCD_H_



#include "ChargeDensity.H"
#include "Functions/FittedFunction.H"

//----------------------------------------------------------------------------------
//
//
class FittedCD
    : public virtual ChargeDensity
    , public virtual FittedFunction
{
public:
    virtual void ReScale      (double factor)=0;  //Ro *= factor
    virtual void ShiftOrigin  (const RVec3& )=0;  //Usefull for single atom charge densities.
    //virtual SMat   GetRepulsion(const BasisSet*) const=0; //unhide
    virtual double GetSelfRepulsion    (                     ) const=0;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual double GetRepulsion        (const FittedFunction*) const=0;  // 1/2 <ro(1) | 1/r12 | ff(2)>

    //Required for creating a polarized CD from and un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
    static  FittedCD*  Factory(std::istream&)        ;
};

#endif //_ChargeDensity_H_
