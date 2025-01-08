// File: RestMass.C  Reast mass c^2 term for the Dirac hamiltonian.

#include "Imp/Hamiltonian/RestMass.H"
#include <BasisSet.H>
#include <TotalEnergy.H>

RestMass::RestMass()
    : HamiltonianTermImp()
{};


HamiltonianTerm::SMat RestMass::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    return bs->GetRestMass();
}

void RestMass::GetEnergy(TotalEnergy& te) const
{
    te.RestMass=CalculateEnergy();
}

std::ostream& RestMass::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "   Rest mass (beta-alpha)*c^2" << std::endl;
    return os;
}


