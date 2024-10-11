// File: ChargeDensity.H  Interface for the charge density category.
#ifndef _ChargeDensity_H_
#define _ChargeDensity_H_

#include "Functions/ScalarFunction.H"
#include "Functions/FittedFunctionClient.H"
#include "Misc/pmstream.h"


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
//Derives from ScalarFunction<double> via PlottableScalarFunction, i.e. it's real.//
//
class ChargeDensity
    : public virtual PMStreamableObject
    , public virtual PlottableScalarFunction
    , public virtual DensityFFClient
{
public:
    typedef SMatrix<double> SMat;
    virtual SMat   GetOverlap  (const IrrepBasisSet*) const=0;
    virtual SMat   GetRepulsion(const IrrepBasisSet*) const=0;
    virtual SMat   GetExchange (const IrrepBasisSet*) const=0;

    virtual double GetEnergy(const HamiltonianTerm*) const=0; //Amounts to Integral(ro*V*d3r);

    virtual double GetTotalCharge      (                     ) const=0;  // <ro>

    virtual bool   IsPolarized  (                           ) const  ;
    virtual void   ReScale      (double factor              )      =0;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      =0;  //Usefull for single atom charge densities.
    virtual void   MixIn        (const ChargeDensity&,double)      =0;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const=0;  //Convergence check.

    virtual double FitGetConstraint(                        ) const;

    virtual std::ostream&  Write(std::ostream&) const=0;
    virtual std::istream&  Read (std::istream&)      =0;
    static  ChargeDensity*  Factory(std::istream&)        ;
};

//---------------------------------------------------------------------------------------
//
//  Store spin up and spin down as a ChargeDensity
//
class PolarizedCD
    : public virtual ChargeDensity
{
public:
    virtual       ChargeDensity* GetChargeDensity(const Spin&)      =0;
    virtual const ChargeDensity* GetChargeDensity(const Spin&) const=0;

    virtual SMat   GetOverlap  (const IrrepBasisSet*) const;
    virtual SMat   GetRepulsion(const IrrepBasisSet*) const;
    virtual SMat   GetExchange (const IrrepBasisSet*) const;

    virtual double GetEnergy(const HamiltonianTerm*) const;

    virtual double GetTotalCharge      (                     ) const;  // <ro>
    virtual double GetTotalSpin        (                     ) const;  // <up>-<down>

    virtual Vector<double> GetRepulsions(const IrrepBasisSet* theFitBasisSet) const;

    virtual bool   IsPolarized  (                           ) const;
    virtual void   ReScale      (double factor              )      ;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;  //Usefull for single atom charge densities.
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //Convergence check.

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;
    virtual void   Eval(const Mesh&, RVec&) const;
};

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
    
    virtual double GetSelfRepulsion    (                     ) const=0;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual double GetRepulsion        (const FittedFunction*) const=0;  // 1/2 <ro(1) | 1/r12 | ff(2)>

    //Required for creating a polarized CD from and un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
    static  FittedCD*  Factory(std::istream&)        ;
};

#include "FunctionsImp/DumpScalarData.H"

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//
class FittedPolarizedCD
    : public virtual FittedCD
    , public virtual PolarizedCD
    , private        DumpScalarData
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

    virtual Vector<double> GetRepulsions(const IrrepBasisSet* theFitBasisSet) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual FittedCD* Clone(        ) const;

private:
    FittedCD* itsSpinUpCD;
    FittedCD* itsSpinDownCD;
};


#endif //_ChargeDensity_H_
