// File: Kinetic.C  Kinetic energy term for the hamiltonian.
module;
#include <iostream>
#include <memory>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Blaze;

namespace qchem::Hamiltonian
{


// The NON-RELATIVISTIC kinetic ENERGY term \f$T=-\tfrac12\nabla^2\f$.  This is the boundary where
// "kinetic energy" is actually born: bs->Kinetic() returns the building block \f$\langle p^2\rangle=
// \langle-\nabla^2\rangle\f$ (no 1/2 -- see BasisSet/Orbital_1E_IBS.C), and the 1/2 is applied HERE.
rsmat_t Kinetic::CalculateMatrix(const robs_t* bs,const Spin&) const
{
    return 0.5*bs->Kinetic();   // T = 1/2 * <p^2>
}

void Kinetic::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    te.Kinetic=cd->DM_Contract(this);
}

std::ostream& Kinetic::Write(std::ostream& os) const
{
    os << "    Kinetic energy Grad^2(r_i)" << std::endl;
    return os;
}


} //namespace
