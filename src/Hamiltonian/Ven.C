// File: Ven.C  Electron-Nuclear potential.



#include "Imp/Hamiltonian/Ven.H"
#include <BasisSet.H>
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include "oml/smatrix.h"
#include "oml/vector.h"
#include "oml/vector3d.h"

Ven::Ven() : HamiltonianTermImp() , theCluster() {};

Ven::Ven(cl_t& cl)
    : HamiltonianTermImp()
    , theCluster(cl)
{
    assert(cl->GetNumAtoms()>0);
};


HamiltonianTerm::SMat Ven::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    return bs->GetNuclear(&*theCluster);
}

void Ven::GetEnergy(TotalEnergy& te) const
{
    te.Een=CalculateEnergy();
}

std::ostream& Ven::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Nuclear-electron potential Zi/|Ri-r| with " << theCluster->GetNumAtoms() << " atoms." << std::endl;
    else
        os << *theCluster;
    return os;
}

std::istream& Ven::Read (std::istream& is)
{
    Cluster* cl=Cluster::Factory(is);
    is >> *cl;
    theCluster.reset(cl);
    return is;
}

