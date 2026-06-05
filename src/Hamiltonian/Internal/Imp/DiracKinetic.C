// File: DiracKinetic.C  Kinetic energy term for the Dirac hamiltonian.
module;
#include <iostream>
#include "blaze/Math.h" 

module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Constants;

namespace qchem::Hamiltonian
{

rsmat_t DiracKinetic::CalculateMatrix(const obs_t* bs,const Spin&) const
{
    // std::cout << "K_dirac=" << bs->Grad2() << std::endl;
    return bs->Kinetic();
}

void DiracKinetic::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Kinetic=cd->DM_Contract(this);
}

std::ostream& DiracKinetic::Write(std::ostream& os) const
{
    os << "    Dirac kinetic energy c*sigma*p" << std::endl;
    return os;
}

} //namespace

