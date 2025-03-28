// File: Kinetic.C  Kinetic energy term for the hamiltonian.

#include "Imp/Hamiltonian/Kinetic.H"
#include <Irrep_BS.H>
#include <TotalEnergy.H>

Kinetic::Kinetic()
    : HamiltonianTermImp()
{};


HamiltonianTerm::SMat Kinetic::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin&) const
{
    return 0.5*bs->Grad2();
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


