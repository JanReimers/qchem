// File: Vxc.C  Hartree-Fock exchange potential

#include "Imp/Hamiltonian/Vxc.H"
#include <HF_IBS.H>
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
    auto hf_bs = dynamic_cast<const TOrbital_HF_IBS<double>*>(bs);
    assert(hf_bs);
    SMat Kab=itsExactCD->GetExchange(hf_bs);
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

