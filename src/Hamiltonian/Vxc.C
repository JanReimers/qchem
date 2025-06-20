// File: Vxc.C  Hartree-Fock exchange potential

#include "Vxc.H"
#include <BasisSet/HF_IBS.H>
#include <ChargeDensity/ChargeDensity.H>
#include <Hamiltonian/TotalEnergy.H>
#include "oml/smatrix.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

Vxc::Vxc() {};

//########################################################################
//
//  Let the charge density do the work.
//

Static_HT::SMat Vxc::CalcMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    auto hf_bs = dynamic_cast<const TOrbital_HF_IBS<double>*>(bs);
    assert(hf_bs);
    SMat Kab=cd->GetExchange(hf_bs);
    return Kab*-0.5;
}
void Vxc::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    te.Exc+=0.5*cd->DM_Contract(this,cd);
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Hartee-Fock exchange potential phi(r_1)*phi(r_2)/r_12" << std::endl;
    return os;
}

