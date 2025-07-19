// File: CDFittedVee.C  Exact Coulomb potential

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

#include "FittedVee.H"
#include <ChargeDensity/Factory.H>
#include <Hamiltonian/TotalEnergy.H>
import qchem.DFT_IBS;
import qchem.ChargeDensity;
import qchem.FittedCD;

FittedVee::FittedVee()
  
{};

FittedVee::FittedVee(bs_t& chargeDensityFitBasisSet, mesh_t&  m, double numElectrons)
    : itsFittedChargeDensity(FittedCD_Factory(chargeDensityFitBasisSet,m,numElectrons))
{
    assert(itsFittedChargeDensity);
};

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro_fit(r2)/r12 .
//              /
//  Where ro is the fitted charge density.
//

Static_HT::SMat FittedVee::CalcMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    auto dft_bs=dynamic_cast<const TOrbital_DFT_IBS<double>*>(bs);
    return itsFittedChargeDensity->GetRepulsion(dft_bs);
}

void FittedVee::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    assert(itsFittedChargeDensity);
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    te.EeeFit    = 0.5*cd->DM_Contract(this,cd);
    te.EeeFitFit = itsFittedChargeDensity->GetSelfRepulsion();
    te.Eee = 2*te.EeeFit - te.EeeFitFit;
}

