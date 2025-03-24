// File: FittedPolarizedCD.C Implementation

#include "Imp/ChargeDensity/FittedPolarizedCD.H"
#include <Spin.H>
#include "oml/vector.h"
#include "oml/smatrix.h"

#include <iostream>
#include <cassert>
//---------------------------------------------------------------------------------
//
//  Construction zone.
//
FittedPolarizedCD::FittedPolarizedCD()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{}; // No UT coverage

FittedPolarizedCD::FittedPolarizedCD(const FittedCD* fcd, double Stotal)
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{
    // No UT coverage
    itsSpinUpCD  =fcd->Clone();
    itsSpinDownCD=fcd->Clone();

    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    double totalCharge=fcd->GetTotalCharge();
    itsSpinUpCD  ->ReScale((totalCharge+Stotal)/(2*totalCharge));
    itsSpinDownCD->ReScale((totalCharge-Stotal)/(2*totalCharge));
};

FittedPolarizedCD::FittedPolarizedCD(FittedCD* up, FittedCD* down)
    : itsSpinUpCD  (up  )
    , itsSpinDownCD(down)
{
    // No UT coverage
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

FittedPolarizedCD::FittedPolarizedCD(const FittedPolarizedCD& pcd)
    : itsSpinUpCD  (pcd.itsSpinUpCD  ->Clone())
    , itsSpinDownCD(pcd.itsSpinDownCD->Clone())
{
    // No UT coverage
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

FittedPolarizedCD::~FittedPolarizedCD()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}


//-------------------------------------------------------------------
//
//  Access to individual components.
//
FittedCD* FittedPolarizedCD::GetChargeDensity(const Spin& S)
{
    // No UT coverage
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    FittedCD* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const FittedCD* FittedPolarizedCD::GetChargeDensity(const Spin& S) const
{
    // No UT coverage
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const FittedCD* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

double FittedPolarizedCD::GetSelfRepulsion() const
{
    // No UT coverage
    return GetRepulsion(this);
}

// <ro(1) | 1/r12 | ff(2)>
double FittedPolarizedCD::GetRepulsion(const FittedFunction* ff) const
{
    // No UT coverage
    double Jup  =itsSpinUpCD  ->GetRepulsion(ff);
    double Jdown=itsSpinDownCD->GetRepulsion(ff);

    return Jup + Jdown;
}

//-------------------------------------------------------------------------
//
//  Fitted function stuff.
//
double FittedPolarizedCD::DoFit(const DensityFFClient& ffc)
{
    // No UT coverage
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const Polarized_CD* polcd=dynamic_cast<const Polarized_CD*>(&ffc);
    assert(polcd);
    double lam_bar=0;
    lam_bar += itsSpinUpCD  ->DoFit(*polcd->GetChargeDensity(Spin::Up  ));
    lam_bar += itsSpinDownCD->DoFit(*polcd->GetChargeDensity(Spin::Down));
    return lam_bar/2.0;
}

double FittedPolarizedCD::DoFit(const ScalarFFClient& ffc)
{
    // No UT coverage
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    double lam_bar=0;
    lam_bar += itsSpinUpCD  ->DoFit(ffc);
    lam_bar += itsSpinDownCD->DoFit(ffc);
    return lam_bar/2.0;
}

Vector<double> FittedPolarizedCD::GetRepulsion3C(const Fit_IBS* fbs) const
{
    // No UT coverage
    return GetChargeDensity(Spin::Up  )->GetRepulsion3C(fbs)
        +  GetChargeDensity(Spin::Down)->GetRepulsion3C(fbs);
    
}

ChargeDensity::SMat FittedPolarizedCD::GetRepulsion(const TOrbital_DFT_IBS<double>* obs) const
{
    // No UT coverage
    return    GetChargeDensity(Spin::Up  )->GetRepulsion(obs)
            + GetChargeDensity(Spin::Down)->GetRepulsion(obs);
}



void FittedPolarizedCD::ReScale(double factor) //Fit *= factor
{
    // No UT coverage
    itsSpinUpCD   -> ReScale(factor);
    itsSpinDownCD -> ReScale(factor);
}

void FittedPolarizedCD::ShiftOrigin(const RVec3& newCenter)
{
    // No UT coverage
    itsSpinUpCD   -> ShiftOrigin(newCenter);
    itsSpinDownCD -> ShiftOrigin(newCenter);
}

void FittedPolarizedCD::FitMixIn(const FittedFunction& ff,double c) // this = this*(1-c) + that*c.
{
    // No UT coverage
    const FittedPolarizedCD* polcd=dynamic_cast<const FittedPolarizedCD*>(&ff);
    assert(polcd);
    itsSpinUpCD   -> FitMixIn(*polcd->itsSpinUpCD  ,c);
    itsSpinDownCD -> FitMixIn(*polcd->itsSpinDownCD,c);
}

double FittedPolarizedCD::FitGetChangeFrom(const FittedFunction& ff) const
{
    // No UT coverage
    return itsSpinUpCD   -> FitGetChangeFrom(ff)
           + itsSpinDownCD -> FitGetChangeFrom(ff);
}

FittedCD* FittedPolarizedCD::Clone() const
{
    return new FittedPolarizedCD(*this);
}

