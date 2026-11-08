// File: Vee_DFT.C  The two DFT Hartree terms, Dunlap and normal.
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

Vee_DFT_Dunlap::Vee_DFT_Dunlap(fbs_t& chargeDensityFitBasisSet, double numElectrons)
{
    // The CD fit basis arrives as the narrow Coulomb-metric (rFIT_CD_ABS) face -- exactly what the
    // density-fitting machinery takes; thread it straight through (no down-cast to the concrete Fit_IBS).
    itsFittedChargeDensity=ChargeDensity::FittedCD_Factory(chargeDensityFitBasisSet,numElectrons);
    assert(itsFittedChargeDensity);
};

Vee_DFT_Dunlap::~Vee_DFT_Dunlap() = default;   // FittedCD is complete here, so the unique_ptr deletes it correctly

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro_fit(r2)/r12 .
//              /
//  Where ro is the fitted charge density.
//

rsmat_t Vee_DFT_Dunlap::MakeMatrix(const robs_t* bs,const Spin& s,DensityID id) const
{
    RequireFitFor(id);   // a READ (V1.41 stage 2): the refresh fitted this density
    auto dft_bs=dynamic_cast<const odftbs_t*>(bs);
    assert(dft_bs);
    return itsFittedChargeDensity->GetRepulsion(dft_bs);
}

// THE EAGER PHASE (R1.0h): fit ONCE, before the block loop, instead of on whichever block asked first.
// Same guard as MakeMatrix's -- newCD is the density-serial test, so driving it here simply moves WHEN the
// one fit happens; a caller outside the phase still gets a correct (lazily fitted) answer.
void Vee_DFT_Dunlap::RefreshForDensity(const rChargeDensity* cd) const
{
    if (cd && newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
}

void Vee_DFT_Dunlap::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
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

Vee_DFT::Vee_DFT(std::unique_ptr<Fitting::GridDensityFit> fit, double numElectrons)
    : itsFit(std::move(fit))
    , itsNumElectrons(numElectrons)
{
    assert(itsFit);
    assert(itsNumElectrons>0);
}

Vee_DFT::~Vee_DFT() = default;

// Every finite density presents its own ProjectedDensity_AO: a matrix its forward, a seed its field.  The fitter asks which.
void Vee_DFT::Fit(const rChargeDensity& cd) const
{
    itsFit->DoFit(dynamic_cast<const Fitting::ProjectedDensity_AO&>(cd));
}

rsmat_t Vee_DFT::MakeMatrix(const robs_t* bs,const Spin&,DensityID id) const
{
    RequireFitFor(id);   // a READ (V1.41 stage 2): the refresh fitted this density
    return itsFit->Repulsion(bs);
}

void Vee_DFT::RefreshForDensity(const rChargeDensity* cd) const
{
    if (cd && newCD(cd)) Fit(*cd);
}

void Vee_DFT::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    if (newCD(cd)) Fit(*cd);
    const double e=itsFit->HartreeEnergy(cd->DM_Contract(this,cd)), dq=itsFit->Integral()-itsNumElectrons;
    te.Add("Eee", e, EnergyRole::Potential);
    te.AddDiagnostic("HartreeChargeError",  dq);
    te.AddDiagnostic("HartreeChargeEnergy", 2*e*dq/itsNumElectrons);
}

} //namespace
