// File: Kinetic.C  Kinetic energy term for the hamiltonian.
module;
#include <iostream>
#include <memory>
#include <vector>
#include "blaze/Math.h" 

module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.IrrepBasisSet;



rsmat_t Kinetic::CalculateMatrix(const ibs_t* bs,const Spin&) const
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


