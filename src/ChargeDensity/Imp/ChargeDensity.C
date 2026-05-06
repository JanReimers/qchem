// File: ChargeDensity.C  Interface for the charge density category.
module;
#include <iostream>
#include <cassert>
#include <blaze/Math.h>
module qchem.ChargeDensity;

import qchem.Symmetry.Spin;

namespace qchem::ChargeDensity
{

//----------------------------------------------------------------------------
//
//  Various integrals.
//
void Polarized_CD::AccumulateDirect(rsmat_t& Jab,const Orbital_HF_IBS<double>* bs) const
{
    GetChargeDensity(Spin::Up  )->AccumulateDirect(Jab,bs);
    GetChargeDensity(Spin::Down)->AccumulateDirect(Jab,bs);
}

void Polarized_CD::AccumulateExchange(rsmat_t& Kab,const Orbital_HF_IBS<double>* bs) const
{
    // No UT coverage
    GetChargeDensity(Spin::Up  )->AccumulateExchange(Kab,bs);
    GetChargeDensity(Spin::Down)->AccumulateExchange(Kab,bs);
}

double Polarized_CD::DM_Contract(const Static_CC* v) const
{
    return GetChargeDensity(Spin::Up  )->DM_Contract(v)+GetChargeDensity(Spin::Down)->DM_Contract(v);
}

double Polarized_CD::DM_Contract(const Dynamic_CC* v,const DM_CD* cd) const
{
    return GetChargeDensity(Spin::Up  )->DM_Contract(v,cd)+GetChargeDensity(Spin::Down)->DM_Contract(v,cd);
}


double Polarized_CD::GetTotalCharge() const
{
    return GetChargeDensity(Spin::Up)->GetTotalCharge() + GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

double Polarized_CD::GetTotalSpin() const
{
    // No UT coverage
    return GetChargeDensity(Spin::Up)->GetTotalCharge() - GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

rvec_t Polarized_CD::GetRepulsion3C(const Fit_IBS* fbs) const
{
    return GetChargeDensity(Spin::Up  )->GetRepulsion3C(fbs)
        +  GetChargeDensity(Spin::Down)->GetRepulsion3C(fbs);
}

//-----------------------------------------------------------------------
//
//  Convergence.
//
void Polarized_CD::MixIn(const DM_CD& cd,double c)
{
    const Polarized_CD* pcd = dynamic_cast<const Polarized_CD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::MixIn could not cast cd" << std::endl;
        exit(-1);
    }
    GetChargeDensity(Spin::Up)  -> MixIn(*pcd->GetChargeDensity(Spin::Up  ),c);
    GetChargeDensity(Spin::Down)-> MixIn(*pcd->GetChargeDensity(Spin::Down),c);
}

double Polarized_CD::GetChangeFrom(const DM_CD& cd) const
{
    const Polarized_CD* pcd = dynamic_cast<const Polarized_CD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::GetChangeFrom could not cast cd" << std::endl;
        exit(-1);
    }
    return GetChargeDensity(Spin::Up)  ->GetChangeFrom(*pcd->GetChargeDensity(Spin::Up  ))
           + GetChargeDensity(Spin::Down)->GetChangeFrom(*pcd->GetChargeDensity(Spin::Down)) ;
}

void Polarized_CD::ReScale(double factor)
{
    // No UT coverage
    GetChargeDensity(Spin::Up)  ->ReScale(factor);
    GetChargeDensity(Spin::Down)->ReScale(factor);
}

//----------------------------------------------------------------------------------
//
//  Real space function stuff.
//
double Polarized_CD::operator()(const rvec3_t& r) const
{
    // No UT coverage
    return (*GetChargeDensity(Spin::Up))(r) + (*GetChargeDensity(Spin::Down))(r);
}

rvec3_t Polarized_CD::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    return GetChargeDensity(Spin::Up)->Gradient(r) + GetChargeDensity(Spin::Down)->Gradient(r);
}

SpinDensity::SpinDensity(DM_CD* up,DM_CD* down)
: itsSpinUpCD  (up  )
, itsSpinDownCD(down)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

SpinDensity::~SpinDensity()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}

double SpinDensity::operator()(const rvec3_t& r) const
{
    // No UT coverage
    return (*itsSpinUpCD)(r) - (*itsSpinDownCD)(r);
}

rvec3_t SpinDensity::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    return itsSpinUpCD->Gradient(r) - itsSpinDownCD->Gradient(r);
}

} //namespace