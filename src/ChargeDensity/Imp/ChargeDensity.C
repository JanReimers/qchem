// File: ChargeDensity.C  Interface for the charge density category.
module;
#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <string>
module qchem.ChargeDensity;
import qchem.Symmetry.Spin;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO (each spin block's AO face)
import qchem.Blaze;

namespace qchem::ChargeDensity
{

//----------------------------------------------------------------------------
//
//  Various integrals.
//

// Coulomb sees the TOTAL density: both spin channels scatter into the same per-irrep Fock blocks.
void Polarized_CD::AccumulateDirectAll(std::vector<rsmat_t>& Jall) const
{
    GetChargeDensity(Spin::Up  )->AccumulateDirectAll(Jall);
    GetChargeDensity(Spin::Down)->AccumulateDirectAll(Jall);
}

// The RHF (unpolarized) exchange term sums K[D_up]+K[D_down] into the same blocks (= K[D_total], then the
// term scales by -1/2).  The polarized term instead drives AccumulateExchangeAll on ONE spin's composite.
void Polarized_CD::AccumulateExchangeAll(std::vector<rsmat_t>& Kall) const
{
    GetChargeDensity(Spin::Up  )->AccumulateExchangeAll(Kall);
    GetChargeDensity(Spin::Down)->AccumulateExchangeAll(Kall);
}

double Polarized_CD::DM_Contract(const Static_CC* v) const
{
    return GetChargeDensity(Spin::Up  )->DM_Contract(v)+GetChargeDensity(Spin::Down)->DM_Contract(v);
}

double Polarized_CD::DM_Contract(const Dynamic_CC* v,const DM_CD* cd) const
{
    return GetChargeDensity(Spin::Up  )->DM_Contract(v,cd)+GetChargeDensity(Spin::Down)->DM_Contract(v,cd);
}

// Coulomb / RHF-exchange energy: the same (total-density) blocks contract with BOTH spin channels
// (= Tr(D_total.B)).  The polarized exchange term instead calls DM_ContractBlocks on ONE spin's composite.
double Polarized_CD::DM_ContractBlocks(const std::map<std::string,rsmat_t>& blocks) const
{
    return GetChargeDensity(Spin::Up  )->DM_ContractBlocks(blocks)+GetChargeDensity(Spin::Down)->DM_ContractBlocks(blocks);
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

rvec_t Polarized_CD::GetRepulsion3C(const BasisSet::FIT_CD_ABS* fbs) const
{
    // The spin blocks are finite (molecular) densities, hence ProjectedDensity_AO -- cross-cast to their AO
    // face and sum the projections (the AO face is no longer a forced base of tDM_CD).
    auto* up=dynamic_cast<const Fitting::ProjectedDensity_AO*>(GetChargeDensity(Spin::Up  ));
    auto* dn=dynamic_cast<const Fitting::ProjectedDensity_AO*>(GetChargeDensity(Spin::Down));
    assert(up && dn && "Polarized_CD spin block is not a ProjectedDensity_AO (finite path)");
    return up->GetRepulsion3C(fbs) + dn->GetRepulsion3C(fbs);
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
    AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
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
    AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
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