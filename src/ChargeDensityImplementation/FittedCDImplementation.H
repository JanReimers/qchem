// File: FittedCDImplementation.H  Fitted charge density.
#ifndef _FittedCDImplementation_H_
#define _FittedCDImplementation_H_



#include "ChargeDensity.H"
#include "FunctionsImp/IntegralConstrainedFF.H"

//---------------------------------------------------------------------------------
//
//  A charge density implemented by fitting the real charge density to an
//  auxillary basis set.
//
template <class T> class FittedCDImplementation
    : public virtual FittedCD
    , public         IntegralConstrainedFF<T>
{
public:
    FittedCDImplementation();
    FittedCDImplementation(const rc_ptr<IrrepBasisSet>&, const rc_ptr<Mesh>& );
    FittedCDImplementation(const rc_ptr<IrrepBasisSet>&, const rc_ptr<Mesh>& , double totalCharge);

    virtual SMat   GetOverlap  (const IrrepBasisSet*) const;
    virtual SMat   GetRepulsion(const IrrepBasisSet*) const;
    virtual SMat   GetExchange (const IrrepBasisSet*) const; //makes no sense

    virtual double GetEnergy(const HamiltonianTerm*) const;

    virtual double GetSelfRepulsion    (                     ) const;  //Does GetRepulsion(*this);
    virtual double GetRepulsion        (const FittedFunction*) const;  // <ro(1) | 1/r12 | f(2)>
    virtual double GetOverlap          (const FittedFunction*) const;  // <ro(1) | f(1)>
    virtual double GetTotalCharge      (                     ) const;

    virtual Vector<double> GetRepulsions(const IrrepBasisSet* theFitBasisSet) const;

    virtual void   ReScale      (double factor              )      ;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //MaxAbs(delta Fit coeff).
    virtual double DoFit        (const ScalarFFClient&)      ;
    virtual double DoFit        (const DensityFFClient&)      ;

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual FittedCD* Clone(        ) const;

private:
    using FittedFunctionImplementation<T>::FitGetCharge;
    using FittedFunctionImplementation<T>::FitGet2CenterOverlap;
    using FittedFunctionImplementation<T>::FitGet2CenterRepulsion;
    using FittedFunctionImplementation<T>::itsFitCoeff;
    using FittedFunctionImplementation<T>::itsBasisSet;
    
    virtual void Eval(const Mesh&, Vector<double>&) const;

    const ChargeDensity* itsExactRep;
    mutable double       itsTotalCharge;
};

#endif //_FittedCDImplementation_H_
