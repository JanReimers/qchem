// File: Vee.C  Electron-Electron Coulomb potential



#include "Imp/Hamiltonian/Vee.H"
#include <ChargeDensity.H>
#include <TotalEnergy.H>
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

Vee::Vee()
    : HamiltonianTermImp()
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

HamiltonianTerm::SMat Vee::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    assert(itsExactCD);
    return itsExactCD->GetRepulsion(bs);
}

void Vee::GetEnergy(TotalEnergy& te) const
{
    assert(itsExactCD);
    te.Eee=0.5*CalculateEnergy();
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

std::ostream& Vee::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}




