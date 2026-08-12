// File: CDFittedVee.C  Exact Coulomb potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.ChargeDensity.Factory;
import qchem.ChargeDensity;
import qchem.FittedCD;
import qchem.Hamiltonian.Types;

namespace qchem::Hamiltonian
{

FittedVee::FittedVee(fbs_t& chargeDensityFitBasisSet, double numElectrons)
{
    // The CD fit basis arrives as the narrow Coulomb-metric (rFIT_CD_ABS) face -- exactly what the
    // density-fitting machinery takes; thread it straight through (no down-cast to the concrete Fit_IBS).
    itsFittedChargeDensity=ChargeDensity::FittedCD_Factory(chargeDensityFitBasisSet,numElectrons);
    assert(itsFittedChargeDensity);
};

FittedVee::~FittedVee() = default;   // FittedCD is complete here, so the unique_ptr deletes it correctly

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro_fit(r2)/r12 .
//              /
//  Where ro is the fitted charge density.
//

rsmat_t FittedVee::MakeMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    auto dft_bs=dynamic_cast<const odftbs_t*>(bs);
    assert(dft_bs);
    return itsFittedChargeDensity->GetRepulsion(dft_bs);
}

void FittedVee::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    assert(itsFittedChargeDensity);
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    // Accumulate through locals: te.Eee is the Dunlap combination of THIS term's two fit pieces,
    // so it must not read the (possibly already accumulated) te.EeeFit/te.EeeFitFit fields.
    double eeeFit   =0.5*cd->DM_Contract(this,cd);
    double eeeFitFit=itsFittedChargeDensity->GetSelfRepulsion();
    te.EeeFit    += eeeFit;
    te.EeeFitFit += eeeFitFit;
    te.Eee       += 2*eeeFit - eeeFitFit;
}

} //namespace
