// File: HartreeFockVxc.C  HartreeFock Coulomb potential



#include "HamiltonianImplementation/HartreeFockVxc.H"
#include "BasisSet/TBasisSet.H"
#include "ChargeDensity/ChargeDensity.H"
#include "Hamiltonian/TotalEnergy.H"
#include "oml/smatrix.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>

HartreeFockVxc::HartreeFockVxc()
    : HamiltonianTermImplementation( )
{};

//########################################################################
//
//  Let the charge density do the work.
//

HamiltonianTerm::SMat HartreeFockVxc::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
{
    assert(itsExactCD);
    SMat Kab=itsExactCD->GetExchange(bs);
    return Kab*-0.5;
}
void HartreeFockVxc::GetEnergy(TotalEnergy& te) const
{
    assert(itsExactCD);
    te.Exc+=0.5*CalculateEnergy();
}

