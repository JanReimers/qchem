// File: ChargeDensity.H  Interface for the charge density category.
#ifndef _ChargeDensity_H_
#define _ChargeDensity_H_

#include <ScalarFunction.H>
#include <FittedFunctionClient.H>


//----------------------------------------------------------------------------------
//
//  Charge density has a simple mandate:
//    1) Provide numerical evluation of ro(r).
//    2) Calculate the Coulomb self energy = sum ni <i(1)|Ro(2)/r12|i(1)> = sum Dab <a(1)|Ro(2)/r12|b(1)>
//    3) Calculate Vcoul(0) = <Ro(r)/r>.
//    4) Calculate the overlap   integrals  < ro(1)| b(1) > for some basis set b.
//    5) Calculate the repulsion integrals  < ro(1)/r12 | b(2) > for some basis set b.
//    6) Calculate the orbital repulsion integrals  < i(1) | ro(2)/r12 | j(1) > for orbitals i,j.
//    7) calculate the self repulsion = 1/2 <ro(1)|1/r12|ro(2)>
//
//
class ChargeDensity
    : public virtual ScalarFunction<double>
    , public virtual DensityFFClient //Fitted function can be fit to this.
{
public:

    virtual double GetEnergy(const HamiltonianTerm*) const=0; //Amounts to Integral(ro*V*d3r);
    virtual double GetTotalCharge      (                     ) const=0;  // <ro>

    virtual void   ReScale      (double factor              )      =0;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      =0;  //Usefull for single atom charge densities.
    virtual void   MixIn        (const ChargeDensity&,double)      =0;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const=0;  //Convergence check.

    virtual double FitGetConstraint(                        ) const;
};



#include <FittedFunction.H>

//----------------------------------------------------------------------------------
//
//
class FittedCD
    : public virtual ChargeDensity
    , public virtual FittedFunction
{
public:
    typedef SMatrix<double> SMat;
    virtual void ReScale      (double factor)=0;  //Ro *= factor
    virtual void ShiftOrigin  (const RVec3& )=0;  //Usefull for single atom charge densities.
    
    virtual double GetSelfRepulsion    (                     ) const=0;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual double GetRepulsion        (const FittedFunction*) const=0;  // 1/2 <ro(1) | 1/r12 | ff(2)>
    virtual SMat   GetRepulsion(const TOrbital_DFT_IBS<double>*) const=0;
    //Required for creating a polarized CD from and un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
};

class Exact_CD : public virtual ChargeDensity
{
public:
    typedef SMatrix<double> SMat;
    virtual SMat   GetRepulsion(const TOrbital_IBS<double>*) const=0;
    virtual SMat   GetExchange (const TOrbital_IBS<double>*) const=0;

};

//---------------------------------------------------------------------------------------
//
//  Store spin up and spin down as a ChargeDensity
//  Generic: Could be fitted or exact.
//
class Polarized_CD
    : public virtual ChargeDensity
{
public:
    virtual       ChargeDensity* GetChargeDensity(const Spin&)      =0;
    virtual const ChargeDensity* GetChargeDensity(const Spin&) const=0;

    virtual double GetEnergy(const HamiltonianTerm*) const;

    virtual double GetTotalCharge      (                     ) const;  // <ro>
    virtual double GetTotalSpin        (                     ) const;  // No UT coverage// <up>-<down>

    virtual Vector<double> GetRepulsion3C(const Fit_IBS*) const;

    virtual void   ReScale      (double factor              )      ;  // No UT coverage//Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;  // No UT coverage//Usefull for single atom charge densities.
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //Convergence check.

    virtual double operator()(const RVec3&) const; // No UT coverage
    virtual RVec3  Gradient  (const RVec3&) const; // No UT coverage
};

class Polarized_Exact_CD 
    : public virtual Polarized_CD
    , public virtual Exact_CD
{
public:
    virtual       Exact_CD* GetChargeDensity(const Spin&)      =0;
    virtual const Exact_CD* GetChargeDensity(const Spin&) const=0;
    virtual SMat   GetRepulsion(const TOrbital_IBS<double>*) const;
    virtual SMat   GetExchange (const TOrbital_IBS<double>*) const; // No UT coverage
};

#endif //_ChargeDensity_H_
