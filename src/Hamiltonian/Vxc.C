// File: Vxc.C  Hartree-Fock exchange potential

#include "Imp/Hamiltonian/Vxc.H"
#include <BasisSet.H>
#include <ChargeDensity.H>
#include <TotalEnergy.H>
#include "oml/smatrix.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

Vxc::Vxc() : HamiltonianTermImp( ) {};

//########################################################################
//
//  Let the charge density do the work.
//

HamiltonianTerm::SMat Vxc::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin&) const
{
    assert(itsExactCD);
    SMat Kab=itsExactCD->GetExchange(bs);
    return Kab*-0.5;
}
void Vxc::GetEnergy(TotalEnergy& te) const
{
    assert(itsExactCD);
    te.Exc+=0.5*CalculateEnergy();
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Hartee-Fock exchange potential phi(r_1)*phi(r_2)/r_12" << std::endl;
    return os;
}

