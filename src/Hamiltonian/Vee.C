// File: Vee.C  Electron-Electron Coulomb potential



#include "Imp/Hamiltonian/Vee.H"
#include <HF_IBS.H>
#include <ChargeDensity.H>
#include <TotalEnergy.H>
#include "oml/smatrix.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

Vee::Vee()
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

Static_HT::SMat Vee::CalculateHamiltonianMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
{
    assert(itsCD);
    auto hf_bs = dynamic_cast<const TOrbital_HF_IBS<double>*>(bs);
    assert(hf_bs);
    return itsCD->GetRepulsion(hf_bs);
}

void Vee::GetEnergy(TotalEnergy& te,const DM_CD* cd) const
{
    assert(itsCD);
    te.Eee=0.5*CalculateEnergy(cd);
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

std::ostream& Vee::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}




