// File: FittedCDImplementation.H  Fitted charge density.
#ifndef _FittedCDImplementation_H_
#define _FittedCDImplementation_H_



#include "ChargeDensity/FittedCD.H"
#include "FunctionsImp/IntegralConstrainedFF.H"
#include "FunctionsImp/DumpScalarData.H"

//---------------------------------------------------------------------------------
//
//  A charge density implemented by fitting the real charge density to an
//  auxillary basis set.
//
template <class T> class FittedCDImplementation
    : public virtual FittedCD
    , public         IntegralConstrainedFF<T>
    , private        DumpScalarData
{
public:
    FittedCDImplementation();
    FittedCDImplementation(const rc_ptr<BasisSet>&);
    FittedCDImplementation(const rc_ptr<BasisSet>&, double totalCharge);

    virtual SMat   GetOverlap  (const BasisSet*) const;
    virtual SMat   GetRepulsion(const BasisSet*) const;
    virtual SMat   GetExchange (const BasisSet*) const; //makes no sense

    virtual double GetEnergy(const HamiltonianTerm*) const;

    virtual double GetSelfRepulsion    (                     ) const;  //Does GetRepulsion(*this);
    virtual double GetRepulsion        (const FittedFunction*) const;  // <ro(1) | 1/r12 | f(2)>
    virtual double GetOverlap          (const FittedFunction*) const;  // <ro(1) | f(1)>
    virtual double GetTotalCharge      (                     ) const;

    virtual void   InjectOverlaps  (FittedFunction*, const BasisSet* theFitBasisSet) const;
    virtual void   InjectRepulsions(FittedFunction*, const BasisSet* theFitBasisSet) const;

    virtual void   ReScale      (double factor              )      ;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //MaxAbs(delta Fit coeff).
    virtual double DoFit        (const FittedFunctionClient&)      ;

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual FittedCD* Clone(        ) const;

private:
    using FittedFunctionImplementation<T>::FitGetCharge;
    using FittedFunctionImplementation<T>::FitGet2CenterOverlap;
    using FittedFunctionImplementation<T>::FitGet2CenterRepulsion;
    virtual void Eval(const Mesh&, Vector<double>&) const;

    const ChargeDensity* itsExactRep;
    mutable double       itsTotalCharge;
};

#endif //_FittedCDImplementation_H_
