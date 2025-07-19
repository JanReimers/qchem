// File: Kinetic.C  Kinetic energy term for the hamiltonian.

#include <iostream>
#include <memory>
#include <vector>

#include "Kinetic.H"
#include <ChargeDensity/ChargeDensity.H>
#include <Hamiltonian/TotalEnergy.H>
import qchem.Irrep_BS;

Kinetic::Kinetic()
    : Static_HT_Imp()
{};


Static_HT::SMat Kinetic::CalculateMatrix(const ibs_t* bs,const Spin&) const
{
    return 0.5*bs->Kinetic();
}

void Kinetic::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Kinetic=cd->DM_Contract(this);
}

std::ostream& Kinetic::Write(std::ostream& os) const
{
    os << "    Kinetic energy Grad^2(r_i)" << std::endl;
    return os;
}


