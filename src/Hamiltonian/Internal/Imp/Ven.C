// File: Ven.C  Electron-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.IrrepBasisSet;

Ven::Ven() : Static_HT_Imp() , theCluster() {};

Ven::Ven(const cl_t& cl)
    : Static_HT_Imp()
    , theCluster(cl)
{
    assert(cl->GetNumAtoms()>0);
};


 SMatrix<double>  Ven::CalculateMatrix(const ibs_t* bs,const Spin&) const
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


