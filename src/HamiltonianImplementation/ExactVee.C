// File: ExactVee.C  Exact Coulomb potential



#include "HamiltonianImplementation/ExactVee.H"
#include "ChargeDensity/ChargeDensity.H"
#include "Hamiltonian/TotalEnergy.H"
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

ExactVee::ExactVee()
    : HamiltonianTermImplementation()
{};


//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro(r2)/r12 .
//              /
//  Where ro is the charge density.
//

HamiltonianTerm::SMat ExactVee::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
{
    assert(itsExactCD);
    return itsExactCD->GetRepulsion(bs);
}

void ExactVee::GetEnergy(TotalEnergy& te) const
{
    assert(itsExactCD);
    te.Eee=0.5*CalculateEnergy();
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

