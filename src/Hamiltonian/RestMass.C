// File: RestMass.C  Reast mass c^2 term for the Dirac hamiltonian.

#include "Imp/Hamiltonian/RestMass.H"
#include <BasisSet.H>
#include <TotalEnergy.H>

RestMass::RestMass()
    : HamiltonianTermImp()
{};


HamiltonianTerm::SMat RestMass::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin&) const
{
    // std::cout << "Rest mass/c^2=" << bs->GetRestMass() << std::endl;
    auto sbs=dynamic_cast<const Orbital_RKB_IBS<double>*>(bs);
    assert(sbs);
    return sbs->RestMass();
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


