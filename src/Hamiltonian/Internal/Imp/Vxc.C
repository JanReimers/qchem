// File: Vxc.C  Hartree-Fock exchange potential
module;
#include <iostream>
#include <vector>
#include <map>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Energy;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  Let the charge density do the work.
//

// Whole-system exchange (doc/ERI4Rework.md §5.4): the density scatters itself across canonical irrep pairs
// (AccumulateExchangeAll -> ScatterBoth on Exchange blocks), so K(j,i) is never built.  The shared
// ContractAll scales every block by Scale() (== itsScale, the Fock K coefficient), so GetMatrix can hand
// back a reference to the already-scaled block.
void Vxc::AccumulateAll(std::vector<rsmat_t>& X,const std::vector<const ohfbs_t*>& hf,const DM_CD* dm) const
{
    dm->AccumulateExchangeAll(X,hf);
}

void Vxc::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    // E_x = 1/2 Tr(D.K_scaled) from this term's own whole-system (already itsScale-scaled) exchange blocks.
    ContractAll(cd);
    te.Exc+=0.5*cd->DM_ContractBlocks(itsJKs);
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    os << "    Hartee-Fock exchange potential phi(r_1)*phi(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
