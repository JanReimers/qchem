// File: PolarizedCD.C  Interface for a spin polarized charge density.


#include "PolarizedCD.H"
#include "FittedPolarizedCD.H"
#include "Misc/Spin.H"
#include "oml/smatrix.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

//----------------------------------------------------------------------------
//
//  Various integrals.
//
ChargeDensity::SMat PolarizedCD::GetOverlap  (const BasisSet* bs) const
{
    return
        GetChargeDensity(Spin::Up  )->GetOverlap(bs) +
        GetChargeDensity(Spin::Down)->GetOverlap(bs);
}

ChargeDensity::SMat PolarizedCD::GetRepulsion(const BasisSet* bs) const
{
//    std::cout.precision(4);
//    std::cout.width(7);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
//    std::cout << "Spin up:" << std::endl;
    SMat Jab_up=GetChargeDensity(Spin::Up  )->GetRepulsion(bs);
//    std::cout << "  Jab_up=:" << Jab_up << std::endl;

//    std::cout << "Spin down:" << std::endl;
    SMat Jab_down=GetChargeDensity(Spin::Down)->GetRepulsion(bs);
//    std::cout << "  Jab_down=:" << Jab_down << std::endl;
//    std::cout << "Jab=:" << Jab_up + Jab_down << std::endl;
    return Jab_up + Jab_down;
}

ChargeDensity::SMat PolarizedCD::GetExchange(const BasisSet* bs) const
{
//    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
//    std::cout << "Spin up:" << std::endl;
    SMat Kab_up=GetChargeDensity(Spin::Up  )->GetExchange(bs);
//    std::cout << "  Kab_up=:" << Kab_up << std::endl;
//
//    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
//    std::cout << "Spin down:" << std::endl;
    SMat Kab_down=GetChargeDensity(Spin::Down)->GetExchange(bs);
//    std::cout << "  Kab_down=:" << Kab_up << std::endl;
    return Kab_up + Kab_down;

//    return
//        GetChargeDensity(Spin::Up  )->GetExchange(bs) +
//        GetChargeDensity(Spin::Down)->GetExchange(bs);
}

double PolarizedCD::GetEnergy(const HamiltonianTerm* v) const
{
    return GetChargeDensity(Spin::Up  )->GetEnergy(v)+GetChargeDensity(Spin::Down)->GetEnergy(v);
}

double PolarizedCD::GetTotalCharge() const
{
    return GetChargeDensity(Spin::Up)->GetTotalCharge() + GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

double PolarizedCD::GetTotalSpin() const
{
    return GetChargeDensity(Spin::Up)->GetTotalCharge() - GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

void PolarizedCD::InjectOverlaps(FittedFunction* ff, const BasisSet* theFitBasisSet) const
{
    FittedPolarizedCD* fpcd=dynamic_cast<FittedPolarizedCD*>(ff);
    if (fpcd)
    {
        FittedCD* up  =dynamic_cast<FittedCD*>(fpcd->GetChargeDensity(Spin::Up  ));
        assert(up);
        FittedCD* down=dynamic_cast<FittedCD*>(fpcd->GetChargeDensity(Spin::Down));
        assert(down);
        GetChargeDensity(Spin::Up  )->InjectOverlaps(up  ,theFitBasisSet);
        GetChargeDensity(Spin::Down)->InjectOverlaps(down,theFitBasisSet);
    }
    else
    {
        GetChargeDensity(Spin::Up  )->InjectOverlaps(ff,theFitBasisSet);
        GetChargeDensity(Spin::Down)->InjectOverlaps(ff,theFitBasisSet);
    }
}

void PolarizedCD::InjectRepulsions(FittedFunction* ff, const BasisSet* theFitBasisSet) const
{
    FittedPolarizedCD* fpcd=dynamic_cast<FittedPolarizedCD*>(ff);
    if (fpcd)
    {
        FittedCD* up  =dynamic_cast<FittedCD*>(fpcd->GetChargeDensity(Spin::Up  ));
        assert(up);
        FittedCD* down=dynamic_cast<FittedCD*>(fpcd->GetChargeDensity(Spin::Down));
        assert(down);
        GetChargeDensity(Spin::Up  )->InjectRepulsions(up  ,theFitBasisSet);
        GetChargeDensity(Spin::Down)->InjectRepulsions(down,theFitBasisSet);
    }
    else
    {
        GetChargeDensity(Spin::Up  )->InjectRepulsions(ff,theFitBasisSet);
        GetChargeDensity(Spin::Down)->InjectRepulsions(ff,theFitBasisSet);
    }
}

//-----------------------------------------------------------------------
//
//  Convergence and origin shifting.
//
void   PolarizedCD::ShiftOrigin(const RVec3& newcenter)
{
    GetChargeDensity(Spin::Up)  ->ShiftOrigin(newcenter) ;
    GetChargeDensity(Spin::Down)->ShiftOrigin(newcenter) ;
}

void PolarizedCD::MixIn(const ChargeDensity& cd,double c)
{
    const PolarizedCD* pcd = dynamic_cast<const PolarizedCD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::MixIn could not cast cd" << std::endl;
        exit(-1);
    }
    GetChargeDensity(Spin::Up)  -> MixIn(*pcd->GetChargeDensity(Spin::Up  ),c);
    GetChargeDensity(Spin::Down)-> MixIn(*pcd->GetChargeDensity(Spin::Down),c);
}

double PolarizedCD::GetChangeFrom(const ChargeDensity& cd) const
{
    const PolarizedCD* pcd = dynamic_cast<const PolarizedCD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::GetChangeFrom could not cast cd" << std::endl;
        exit(-1);
    }
    return GetChargeDensity(Spin::Up)  ->GetChangeFrom(*pcd->GetChargeDensity(Spin::Up  ))
           + GetChargeDensity(Spin::Down)->GetChangeFrom(*pcd->GetChargeDensity(Spin::Down)) ;
}

void PolarizedCD::ReScale(double factor)
{
    GetChargeDensity(Spin::Up)  ->ReScale(factor);
    GetChargeDensity(Spin::Down)->ReScale(factor);
}


bool PolarizedCD::IsPolarized() const
{
    return true;
}

//----------------------------------------------------------------------------------
//
//  Real space function stuff.
//
double PolarizedCD::operator()(const RVec3& r) const
{
    return (*GetChargeDensity(Spin::Up))(r) + (*GetChargeDensity(Spin::Down))(r);
}

ChargeDensity::RVec3 PolarizedCD::Gradient  (const RVec3& r) const
{
    return GetChargeDensity(Spin::Up)->Gradient(r) + GetChargeDensity(Spin::Down)->Gradient(r);
}

void PolarizedCD::Eval(const Mesh& m, Vec& v) const
{
    GetChargeDensity(Spin::Up)  ->Eval(m,v);
    GetChargeDensity(Spin::Down)->Eval(m,v);
}

