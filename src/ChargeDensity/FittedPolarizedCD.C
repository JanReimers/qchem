// File: FittedPolarizedCD.C Implementation

#include "Imp/ChargeDensity/FittedPolarizedCD.H"
#include <Spin.H>
#include "oml/vector.h"
#include <iostream>
#include <cassert>
//---------------------------------------------------------------------------------
//
//  Construction zone.
//
FittedPolarizedCD::FittedPolarizedCD()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{};

FittedPolarizedCD::FittedPolarizedCD(const ChargeDensity* unpolcd, double Stotal)
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{
    const FittedCD* fcd=dynamic_cast<const FittedCD*>(unpolcd);
    if (!fcd)
    {
        std::cerr << "FittedPolarizedCD::FittedPolarizedCD could cast input charge density" << std::endl;
        exit(-1);
    }
    itsSpinUpCD  =fcd->Clone();
    itsSpinDownCD=fcd->Clone();

    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    double totalCharge=unpolcd->GetTotalCharge();
    itsSpinUpCD  ->ReScale((totalCharge+Stotal)/(2*totalCharge));
    itsSpinDownCD->ReScale((totalCharge-Stotal)/(2*totalCharge));
};

FittedPolarizedCD::FittedPolarizedCD(ChargeDensity* up, ChargeDensity* down)
    : itsSpinUpCD  (dynamic_cast<FittedCD*>(up  ))
    , itsSpinDownCD(dynamic_cast<FittedCD*>(down))
{
    if (!itsSpinUpCD)
    {
        std::cerr << "FittedPolarizedCD::FittedPolarizedCD could cast up charge density" << std::endl;
        exit(-1);
    }
    if (!itsSpinDownCD)
    {
        std::cerr << "FittedPolarizedCD::FittedPolarizedCD could cast down charge density" << std::endl;
        exit(-1);
    }
};

FittedPolarizedCD::FittedPolarizedCD(const FittedPolarizedCD& pcd)
    : itsSpinUpCD  (pcd.itsSpinUpCD  ->Clone())
    , itsSpinDownCD(pcd.itsSpinDownCD->Clone())
{
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
ChargeDensity* FittedPolarizedCD::GetChargeDensity(const Spin& S)
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    ChargeDensity* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const ChargeDensity* FittedPolarizedCD::GetChargeDensity(const Spin& S) const
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const ChargeDensity* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

double FittedPolarizedCD::GetSelfRepulsion() const
{
    return GetRepulsion(this);
}

// <ro(1) | 1/r12 | ff(2)>
double FittedPolarizedCD::GetRepulsion(const FittedFunction* ff) const
{
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
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const PolarizedCD* polcd=dynamic_cast<const PolarizedCD*>(&ffc);
    assert(polcd);
    double lam_bar=0;
    lam_bar += itsSpinUpCD  ->DoFit(*polcd->GetChargeDensity(Spin::Up  ));
    lam_bar += itsSpinDownCD->DoFit(*polcd->GetChargeDensity(Spin::Down));
    return lam_bar/2.0;
}

double FittedPolarizedCD::DoFit(const ScalarFFClient& ffc)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    double lam_bar=0;
    lam_bar += itsSpinUpCD  ->DoFit(ffc);
    lam_bar += itsSpinDownCD->DoFit(ffc);
    return lam_bar/2.0;
}

Vector<double> FittedPolarizedCD::GetRepulsion3C(const IrrepBasisSet* theFitBasisSet) const
{
    return GetChargeDensity(Spin::Up  )->GetRepulsion3C(theFitBasisSet)
        +  GetChargeDensity(Spin::Down)->GetRepulsion3C(theFitBasisSet);
    
}


void FittedPolarizedCD::ReScale(double factor) //Fit *= factor
{
    itsSpinUpCD   -> ReScale(factor);
    itsSpinDownCD -> ReScale(factor);
}

void FittedPolarizedCD::ShiftOrigin(const RVec3& newCenter)
{
    itsSpinUpCD   -> ShiftOrigin(newCenter);
    itsSpinDownCD -> ShiftOrigin(newCenter);
}

void FittedPolarizedCD::FitMixIn(const FittedFunction& ff,double c) // this = this*(1-c) + that*c.
{
    const FittedPolarizedCD* polcd=dynamic_cast<const FittedPolarizedCD*>(&ff);
    assert(polcd);
    itsSpinUpCD   -> FitMixIn(*polcd->itsSpinUpCD  ,c);
    itsSpinDownCD -> FitMixIn(*polcd->itsSpinDownCD,c);
}

double FittedPolarizedCD::FitGetChangeFrom(const FittedFunction& ff) const
{
    return itsSpinUpCD   -> FitGetChangeFrom(ff)
           + itsSpinDownCD -> FitGetChangeFrom(ff);
}

//--------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& FittedPolarizedCD::Write(std::ostream& os) const
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    return os << itsSpinUpCD << itsSpinDownCD;
}

std::istream& FittedPolarizedCD::Read (std::istream& is)
{
    delete itsSpinUpCD;
    itsSpinUpCD=FittedCD::Factory(is);
    assert(itsSpinUpCD);
    is >> *itsSpinUpCD;

    delete itsSpinDownCD;
    itsSpinDownCD=FittedCD::Factory(is);
    assert(itsSpinDownCD);
    is >> *itsSpinDownCD;

    return is;
}

FittedCD* FittedPolarizedCD::Clone() const
{
    return new FittedPolarizedCD(*this);
}

