// File: Kinetic.C  Kinetic energy term for the hamiltonian.

#include "Imp/Hamiltonian/Kinetic.H"
#include <BasisSet.H>
#include <TotalEnergy.H>

Kinetic::Kinetic()
    : HamiltonianTermImp()
{};


HamiltonianTerm::SMat Kinetic::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    //std::cout << "K=" << bs->GetKinetic() << std::endl;
    
    return bs->GetKinetic();
}

void Kinetic::GetEnergy(TotalEnergy& te) const
{
    te.Kinetic=CalculateEnergy();
}

std::ostream& Kinetic::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Kinetic energy Grad^2(r_i)" << std::endl;
    return os;
}


