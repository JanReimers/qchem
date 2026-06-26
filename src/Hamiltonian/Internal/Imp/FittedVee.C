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
import qchem.BasisSet.Fit_IBS;   // Fit_IBS (the full fit basis FittedCD holds; recovered from the narrow CD face)

namespace qchem::Hamiltonian
{

FittedVee::FittedVee(bs_t& chargeDensityFitBasisSet, double numElectrons)
{
    // The FittedCD machinery (qcChargeDensity) holds the full Fit_IBS; recover it from the narrow CD face
    // we were handed for type safety.  Always succeeds: every fit basis the creators build IS a Fit_IBS.
    auto fbs=std::dynamic_pointer_cast<const BasisSet::Fit_IBS>(chargeDensityFitBasisSet);
    assert(fbs && "FittedVee: the CD fit basis must be a concrete Fit_IBS");
    itsFittedChargeDensity=ChargeDensity::FittedCD_Factory(fbs,numElectrons);
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

rsmat_t FittedVee::CalcMatrix(const obs_t* bs,const Spin& s,const DM_CD* cd) const
{
    if (newCD(cd)) itsFittedChargeDensity->DoFit(*cd);
    auto dft_bs=dynamic_cast<const odftbs_t*>(bs);
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

} //namespace
