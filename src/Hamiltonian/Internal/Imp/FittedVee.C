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

// THE EAGER PHASE (R1.0h): fit ONCE, before the block loop, instead of on whichever block asked first.
// Same guard as MakeMatrix's -- newCD is the density-serial test, so driving it here simply moves WHEN the
// one fit happens; a caller outside the phase still gets a correct (lazily fitted) answer.
void FittedVee::RefreshForDensity(const rChargeDensity* cd) const
{
    if (cd && newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
}

void FittedVee::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    assert(itsFittedChargeDensity);
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    // The Dunlap combination of THIS term's two fit pieces is the contribution; the pieces themselves are
    // diagnostics (reported, never summed).  Tr(D V) is left ABSENT: under density mixing the fitted V_H's
    // expectation in the Fock is Tr(D_out V[rho_in]) -- the Harris-Foulkes subtlety -- not 2 eeeFit.
    double eeeFit   =0.5*cd->DM_Contract(this,cd);
    double eeeFitFit=itsFittedChargeDensity->GetSelfRepulsion();
    te.AddDiagnostic("EeeFit",    eeeFit);
    te.AddDiagnostic("EeeFitFit", eeeFitFit);
    te.Add("Eee", 2*eeeFit - eeeFitFit, EnergyRole::Potential);
}

} //namespace
