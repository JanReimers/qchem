// File: ChargeDensity.C  Interface for the charge density category.
#include "ChargeDensity.H"

// double ChargeDensity::FitGetConstraint  () const
// {
//     return  GetTotalCharge();
// }

#include "oml/smatrix.h"
#include <Spin.H>
#include "oml/vector.h"
//----------------------------------------------------------------------------
//
//  Various integrals.
//
Exact_CD::SMat Polarized_CD::GetRepulsion(const TOrbital_HF_IBS<double>* bs) const
{
    SMat Jab_up=GetChargeDensity(Spin::Up  )->GetRepulsion(bs);
    SMat Jab_down=GetChargeDensity(Spin::Down)->GetRepulsion(bs);
    return Jab_up + Jab_down;
}

Exact_CD::SMat Polarized_CD::GetExchange(const TOrbital_HF_IBS<double>* bs) const
{
    // No UT coverage
    SMat Kab_up=GetChargeDensity(Spin::Up  )->GetExchange(bs);
    SMat Kab_down=GetChargeDensity(Spin::Down)->GetExchange(bs);
    return Kab_up + Kab_down;
}

double Polarized_CD::GetEnergy(const HamiltonianTerm* v) const
{
    return GetChargeDensity(Spin::Up  )->GetEnergy(v)+GetChargeDensity(Spin::Down)->GetEnergy(v);
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

Vector<double> Polarized_CD::GetRepulsion3C(const Fit_IBS* fbs) const
{
    return GetChargeDensity(Spin::Up  )->GetRepulsion3C(fbs)
        +  GetChargeDensity(Spin::Down)->GetRepulsion3C(fbs);
}

//-----------------------------------------------------------------------
//
//  Convergence and origin shifting.
//
void   Polarized_CD::ShiftOrigin(const RVec3& newcenter)
{
    // No UT coverage
    GetChargeDensity(Spin::Up)  ->ShiftOrigin(newcenter) ;
    GetChargeDensity(Spin::Down)->ShiftOrigin(newcenter) ;
}

void Polarized_CD::MixIn(const Exact_CD& cd,double c)
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

double Polarized_CD::GetChangeFrom(const Exact_CD& cd) const
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
double Polarized_CD::operator()(const RVec3& r) const
{
    // No UT coverage
    return (*GetChargeDensity(Spin::Up))(r) + (*GetChargeDensity(Spin::Down))(r);
}

Exact_CD::RVec3 Polarized_CD::Gradient  (const RVec3& r) const
{
    // No UT coverage
    return GetChargeDensity(Spin::Up)->Gradient(r) + GetChargeDensity(Spin::Down)->Gradient(r);
}


