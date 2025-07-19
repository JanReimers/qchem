// File: Ven.C  Electron-Nuclear potential.


#include <iostream>
#include <cassert>
#include <memory>
#include <vector>

#include "Ven.H"
#include <Hamiltonian/TotalEnergy.H>
#include <ChargeDensity/ChargeDensity.H>
import qchem.Irrep_BS;

Ven::Ven() : Static_HT_Imp() , theCluster() {};

Ven::Ven(const cl_t& cl)
    : Static_HT_Imp()
    , theCluster(cl)
{
    assert(cl->GetNumAtoms()>0);
};


Static_HT::SMat Ven::CalculateMatrix(const ibs_t* bs,const Spin&) const
{
    // std::cout << "Ven=" << bs->Nuclear(&*theCluster) << std::endl;
    return bs->Nuclear(&*theCluster);
}

void Ven::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Een=cd->DM_Contract(this);
}

std::ostream& Ven::Write(std::ostream& os) const
{
    os << "    Nuclear-electron potential Zi/|Ri-r| with " << theCluster->GetNumAtoms() << " atoms." << std::endl;
    return os;
}


