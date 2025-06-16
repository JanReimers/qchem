// File: Ven.C  Electron-Nuclear potential.



#include "Imp/Hamiltonian/Ven.H"
#include <BasisSet/Irrep_BS.H>
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include "oml/smatrix.h"
#include "oml/vector.h"
#include "oml/vector3d.h"

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
    if (StreamableObject::Pretty())
        os << "    Nuclear-electron potential Zi/|Ri-r| with " << theCluster->GetNumAtoms() << " atoms." << std::endl;
    else
        os << *theCluster;
    return os;
}


