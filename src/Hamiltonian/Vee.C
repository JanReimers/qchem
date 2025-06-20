// File: Vee.C  Electron-Electron Coulomb potential



#include "Vee.H"
#include <BasisSet/HF_IBS.H>
#include <ChargeDensity/ChargeDensity.H>
#include <Hamiltonian/TotalEnergy.H>
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

Static_HT::SMat Vee::CalcMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    auto hf_bs = dynamic_cast<const TOrbital_HF_IBS<double>*>(bs);
    assert(hf_bs);
    return cd->GetRepulsion(hf_bs);
}

void Vee::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    te.Eee=0.5*cd->DM_Contract(this,cd);
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

std::ostream& Vee::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}




