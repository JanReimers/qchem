// File: CDFittedVee.C  Exact Coulomb potential



#include "Imp/Hamiltonian/FittedVee.H"
#include <ChargeDensity.H>
#include <DFT_IBS.H>
#include "Imp/ChargeDensity/FittedCD.H"
#include <TotalEnergy.H>
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

FittedVee::FittedVee()
    : HamiltonianTermImp()
{};

FittedVee::FittedVee(bs_t& chargeDensityFitBasisSet, mesh_t&  m, double numElectrons)
    : HamiltonianTermImp()
    , itsFittedChargeDensity(new FittedCDImp<double>(chargeDensityFitBasisSet,m,numElectrons))
{
    assert(itsFittedChargeDensity);
};

void FittedVee::UseChargeDensity(const Exact_CD* cd)
{
    assert(cd);
    HamiltonianTermImp::UseChargeDensity(cd);
    itsFittedChargeDensity->DoFit(*cd);
}
//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro_fit(r2)/r12 .
//              /
//  Where ro is the fitted charge density.
//

HamiltonianTerm::SMat FittedVee::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin&) const
{
    auto dft_bs=dynamic_cast<const TOrbital_DFT_IBS<double>*>(bs);
    return itsFittedChargeDensity->GetRepulsion(dft_bs);
}

void FittedVee::GetEnergy(TotalEnergy& te) const
{
    assert(itsFittedChargeDensity);
    assert(itsExactCD);
    te.EeeFit    = 0.5*CalculateEnergy();
    te.EeeFitFit = itsFittedChargeDensity->GetSelfRepulsion();
    te.Eee = 2*te.EeeFit - te.EeeFitFit;
}

