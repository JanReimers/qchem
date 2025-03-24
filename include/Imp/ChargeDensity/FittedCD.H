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
    FittedCDImp(); // No UT coverage
    FittedCDImp(bs_t&, mesh_t&); // No UT coverage
    FittedCDImp(bs_t&, mesh_t&, double totalCharge);

    virtual SMat   GetRepulsion(const TOrbital_DFT_IBS<double>*) const;
    
    // virtual double GetEnergy(const HamiltonianTerm*) const; // No UT coverage

    virtual double GetSelfRepulsion    (                     ) const;  //Does GetRepulsion(*this);
    
    virtual double DoFit        (const DensityFFClient&)      ;
    virtual double DoFit        (const Exact_CD&);

    virtual double operator()(const RVec3&) const; // No UT coverage
    virtual Vec3   Gradient  (const RVec3&) const; // No UT coverage

    virtual FittedCD* Clone(        ) const;

private:
    using FittedFunctionImp<T>::FitGetCharge;
    using FittedFunctionImp<T>::FitGet2CenterOverlap;
    using FittedFunctionImp<T>::FitGet2CenterRepulsion;
    using FittedFunctionImp<T>::itsFitCoeff;
    using FittedFunctionImp<T>::itsBasisSet;
    
    virtual void Eval(const Mesh&, Vector<double>&) const; // No UT coverage

    const Exact_CD* itsExactRep;
    mutable double       itsTotalCharge;
};

#endif //_FittedCDImp_H_
