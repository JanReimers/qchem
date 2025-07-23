// File: Vee.C  Electron-Electron Coulomb potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.HF_IBS;
import qchem.ChargeDensity;
import qchem.Energy;

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

 SMatrix<double>  Vee::CalcMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
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
    os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}




