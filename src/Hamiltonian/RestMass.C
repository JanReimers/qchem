// File: RestMass.C  Reast mass c^2 term for the Dirac hamiltonian.

#include "Imp/Hamiltonian/RestMass.H"
#include "Imp/Misc/DFTDefines.H"
#include <DHF_IBS.H>
#include <TotalEnergy.H>

RestMass::RestMass()
    : Static_HT_Imp()
{};


Static_HT::SMat RestMass::CalculateHamiltonianMatrix(const ibs_t* bs,const Spin&) const
{
    static const double f=-2.0*c_light*c_light;
    // std::cout << "Rest mass/c^2=" << bs->GetRestMass() << std::endl;
    auto sbs=dynamic_cast<const Orbital_RKB_IBS<double>*>(bs);
    assert(sbs);
    return f*sbs->RestMass();
}

void RestMass::GetEnergy(TotalEnergy& te,const Exact_CD* cd) const
{
    te.RestMass=CalculateEnergy(cd);
}

std::ostream& RestMass::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "   Rest mass (beta-alpha)*c^2" << std::endl;
    return os;
}


