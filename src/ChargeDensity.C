// File: ChargeDensity.C  Interface for the charge density category.
#include "ChargeDensity.H"

double ChargeDensity::FitGetConstraint  () const
{
    return  GetTotalCharge();
}

bool ChargeDensity::IsPolarized() const
{
    return false;
}

#include "oml/smatrix.h"
#include <Spin.H>
#include "oml/vector.h"
//----------------------------------------------------------------------------
//
//  Various integrals.
//
ChargeDensity::SMat PolarizedCD::GetRepulsion(const IrrepBasisSet* bs) const
{
    SMat Jab_up=GetChargeDensity(Spin::Up  )->GetRepulsion(bs);
    SMat Jab_down=GetChargeDensity(Spin::Down)->GetRepulsion(bs);
    return Jab_up + Jab_down;
}

ChargeDensity::SMat PolarizedCD::GetExchange(const IrrepBasisSet* bs) const
{
    SMat Kab_up=GetChargeDensity(Spin::Up  )->GetExchange(bs);
    SMat Kab_down=GetChargeDensity(Spin::Down)->GetExchange(bs);
    return Kab_up + Kab_down;
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

Vector<double> PolarizedCD::GetRepulsions(const IrrepBasisSet* theFitBasisSet) const
{
    return GetChargeDensity(Spin::Up  )->GetRepulsions(theFitBasisSet)
        +  GetChargeDensity(Spin::Down)->GetRepulsions(theFitBasisSet);
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


