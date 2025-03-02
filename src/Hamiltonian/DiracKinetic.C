// File: DiracKinetic.C  Kinetic energy term for the Dirac hamiltonian.

#include "Imp/Hamiltonian/DiracKinetic.H"
#include "Imp/Misc/DFTDefines.H"
#include <BasisSet.H>
#include <TotalEnergy.H>

DiracKinetic::DiracKinetic()
    : HamiltonianTermImp()
{};


HamiltonianTerm::SMat DiracKinetic::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin&) const
{
    // std::cout << "K_dirac/c=" << bs->Kinetic() << std::endl;
    return c_light*bs->Kinetic();
}

void DiracKinetic::GetEnergy(TotalEnergy& te) const
{
    te.Kinetic=0.5*CalculateEnergy();
}

std::ostream& DiracKinetic::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Dirac kinetic energy c*sigma*p" << std::endl;
    return os;
}


