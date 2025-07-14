// File: DiracKinetic.C  Kinetic energy term for the Dirac hamiltonian.

#include "DiracKinetic.H"
#include <BasisSet/Irrep_BS.H>
#include <ChargeDensity/ChargeDensity.H>
#include <Hamiltonian/TotalEnergy.H>

import Common.Constants;

DiracKinetic::DiracKinetic()
    : Static_HT_Imp()
{};


Static_HT::SMat DiracKinetic::CalculateMatrix(const ibs_t* bs,const Spin&) const
{
    // std::cout << "K_dirac/c=" << bs->Grad2() << std::endl;
    return c_light*bs->Kinetic();
}

void DiracKinetic::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Kinetic=0.5*cd->DM_Contract(this);
}

std::ostream& DiracKinetic::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Dirac kinetic energy c*sigma*p" << std::endl;
    return os;
}


