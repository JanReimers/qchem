// File: FittedCDImp.H  Fitted charge density.
#ifndef _FittedCDImp_H_
#define _FittedCDImp_H_



#include <ChargeDensity.H>
#include "Imp/Fitting/IntegralConstrainedFF.H"

//---------------------------------------------------------------------------------
//
//  A charge density implemented by fitting the real charge density to an
//  auxillary basis set.
//
template <class T> class FittedCDImp
    : public virtual FittedCD
    , public         IntegralConstrainedFF<T>
{
     typedef typename IntegralConstrainedFF<T>::mesh_t mesh_t;
     typedef typename IntegralConstrainedFF<T>::bs_t   bs_t;
public:
    FittedCDImp();
    FittedCDImp(bs_t&, mesh_t&);
    FittedCDImp(bs_t&, mesh_t&, double totalCharge);

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
    using FittedFunctionImp<T>::FitGetCharge;
    using FittedFunctionImp<T>::FitGet2CenterOverlap;
    using FittedFunctionImp<T>::FitGet2CenterRepulsion;
    using FittedFunctionImp<T>::itsFitCoeff;
    using FittedFunctionImp<T>::itsBasisSet;
    
    virtual void Eval(const Mesh&, Vector<double>&) const;

    const ChargeDensity* itsExactRep;
    mutable double       itsTotalCharge;
};

#endif //_FittedCDImp_H_
