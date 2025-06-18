// File: ChargeDensity.H  Interface for the charge density category.
#ifndef _ChargeDensity_H_
#define _ChargeDensity_H_

#include <Symmetry/fwd.H>
#include <Mesh/ScalarFunction.H>
#include <ChargeDensity/fwd.H>
#include <Fitting/FittedFunctionClient.H>

//
//
//  These little interfaces allow us to invert a dependency with Hamiltonian Terms.
class Static_CC //Contract client for static Ham terms.
{
public:
    virtual const SMatrix<double>& GetMatrix(const TOrbital_IBS<double>*,const Spin&) const=0;    
};

class Dynamic_CC //Contract client for dynamic (CD dependent) Ham terms.
{
public:
    virtual const SMatrix<double>& GetMatrix(const TOrbital_IBS<double>*,const Spin&,const DM_CD*) const=0;    
};

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
//  This is the interface for a charge density representation based on the density matrix.
//
class DM_CD 
: public virtual ScalarFunction<double>
, public virtual DensityFFClient //Fitted function can be fit to this.
{
public:
    virtual double DM_Contract(const Static_CC*) const=0; //Amounts to Integral(ro*V*d3r);
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const=0; //Amounts to Integral(ro*V(ro)*d3r);

    virtual void   ReScale      (double factor         )      =0;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&          )      =0;  //Usefull for single atom charge densities.
    virtual void   MixIn        (const DM_CD&,double)      =0;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const=0;  //Convergence check.

    virtual double GetTotalCharge  () const=0;  // <ro>
    virtual double FitGetConstraint() const {return  GetTotalCharge();}

    typedef SMatrix<double> SMat;
    virtual SMat   GetRepulsion(const TOrbital_HF_IBS<double>*) const=0;
    virtual SMat   GetExchange (const TOrbital_HF_IBS<double>*) const=0;

};

//---------------------------------------------------------------------------------------
//
//  Store spin up and spin down as a ChargeDensity
//  Generic: Could be fitted or exact.
//
class Polarized_CD
    : public virtual DM_CD
{
public:
    virtual       DM_CD* GetChargeDensity(const Spin&)      =0;
    virtual const DM_CD* GetChargeDensity(const Spin&) const=0;

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;

    virtual double GetTotalCharge() const;  // <ro>
    virtual double GetTotalSpin  () const;  // No UT coverage// <up>-<down>

    virtual Vector<double> GetRepulsion3C(const Fit_IBS*) const;
    virtual SMat   GetRepulsion(const TOrbital_HF_IBS<double>*) const;
    virtual SMat   GetExchange (const TOrbital_HF_IBS<double>*) const; 

    virtual void   ReScale      (double factor              )      ;  // No UT coverage//Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;  // No UT coverage//Usefull for single atom charge densities.
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //Convergence check.

    virtual double operator()(const RVec3&) const; // No UT coverage
    virtual RVec3  Gradient  (const RVec3&) const; // No UT coverage
};

class SpinDensity : public virtual ScalarFunction<double>
{
public:
    SpinDensity(DM_CD* up,DM_CD* down);
    ~SpinDensity();
    virtual double operator()(const RVec3&) const; // No UT coverage
    virtual RVec3  Gradient  (const RVec3&) const; // No UT coverage
private:
    DM_CD* itsSpinUpCD;
    DM_CD* itsSpinDownCD;
};


#endif //_ChargeDensity_H_
