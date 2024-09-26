// File: Kinetic.C  Kinetic energy term for the hamiltonian.

#include "HamiltonianImplementation/Kinetic.H"
#include "BasisSet/BasisSet.H"
#include "Hamiltonian/TotalEnergy.H"

Kinetic::Kinetic()
    : HamiltonianTermImplementation()
{};


HamiltonianTerm::SMat Kinetic::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
{
    return bs->GetKinetic();
}

void Kinetic::GetEnergy(TotalEnergy& te) const
{
    te.Kinetic=CalculateEnergy();
}

