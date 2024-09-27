// File: ChargeDensity.H  Interface for the charge density category.
#ifndef _ChargeDensity_H_
#define _ChargeDensity_H_

#include "Functions/PlottableScalarFunction.H"
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
    , public virtual FittedFunctionClient
{
public:
    typedef SMatrix<double> SMat;
    virtual SMat   GetOverlap  (const BasisSet*) const=0;
    virtual SMat   GetRepulsion(const BasisSet*) const=0;
    virtual SMat   GetExchange (const BasisSet*) const=0;

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

#endif //_ChargeDensity_H_
