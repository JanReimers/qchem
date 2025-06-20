// File: FittedCDImp.H  Fitted charge density.
#ifndef _FittedCDImp_H_
#define _FittedCDImp_H_



#include <ChargeDensity/FittedCD.H>
#include "Imp/Fitting/IntegralConstrainedFF.H"

//---------------------------------------------------------------------------------
//
//  A charge density implemented by fitting the real charge density to an
//  auxillary basis set.
//
template <class T> class FittedCDImp
    : public virtual FittedCD
    , public         IntegralConstrainedFF<T> //Pick up the DoFit method.
{
     typedef typename IntegralConstrainedFF<T>::mesh_t mesh_t;
     typedef typename IntegralConstrainedFF<T>::bs_t   bs_t;
public:
    FittedCDImp(bs_t&, mesh_t&, double totalCharge);

    virtual SMat   GetRepulsion    (const TOrbital_DFT_IBS<double>*) const;
    virtual double GetSelfRepulsion(                      ) const;  //Does GetRepulsion(*this);
    
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
};

#endif //_FittedCDImp_H_
