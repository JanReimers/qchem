// File: CDFittedVee.C  Exact Coulomb potential



#include "Imp/Hamiltonian/FittedVee.H"
#include "ChargeDensity.H"
#include "ChargeDensityImplementation/FittedCDImplementation.H"
#include "TotalEnergy.H"
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

FittedVee::FittedVee()
    : HamiltonianTermImp()
{};

FittedVee::FittedVee(const rc_ptr<IrrepBasisSet>& chargeDensityFitBasisSet, const rc_ptr<Mesh>&  m, double numElectrons)
    : HamiltonianTermImp()
    , itsFittedChargeDensity(new FittedCDImplementation<double>(chargeDensityFitBasisSet,m,numElectrons))
{
    assert(itsFittedChargeDensity);
};

void FittedVee::UseChargeDensity(const ChargeDensity* cd)
{
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

HamiltonianTerm::SMat FittedVee::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    const ChargeDensity* cd=itsFittedChargeDensity;
    return cd->GetRepulsion(bs);
}

void FittedVee::GetEnergy(TotalEnergy& te) const
{
    assert(itsFittedChargeDensity);
    assert(itsExactCD);
    te.EeeFit    = 0.5*CalculateEnergy();
    te.EeeFitFit = itsFittedChargeDensity->GetSelfRepulsion();
    te.Eee = 2*te.EeeFit - te.EeeFitFit;
}

