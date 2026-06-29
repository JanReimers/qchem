// File: Vee.C  Electron-Electron Coulomb potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro(r2)/r12 .
//              /
//  Where ro is the charge density.
//

rsmat_t Vee::CalcMatrix(const obs_t* bs,const Spin&,const rChargeDensity* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    auto hf_bs = dynamic_cast<const ohfbs_t*>(bs);
    assert(hf_bs);
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);   // HF J needs the density matrix (not a fit seed)
    assert(dm && "Vee (HF Coulomb): density must be a DM_CD");
    rsmat_t Jab=blazem::zero<double>(bs->GetNumFunctions());
    dm->AccumulateDirect(Jab,hf_bs);
    return Jab;
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

} //namespace
