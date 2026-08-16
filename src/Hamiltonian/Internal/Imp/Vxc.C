// File: Vxc.C  Hartree-Fock exchange potential
module;
#include <iostream>
#include <stdexcept>
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
void Vxc::AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const
{
    // V1.6: ask for the WHOLE-SYSTEM exact-exchange face rather than calling through the general density
    // face.  Before, a density that could not span the irreps (a bare leaf) hit an assert-only default --
    // a silent NO-OP under -DNDEBUG, i.e. a zeroed E and a wrong Fock in the build we ship.
    auto* sys=dynamic_cast<const qchem::ChargeDensity::tHF_System_CD<double>*>(dm);
    if (!sys) throw std::runtime_error("HF term: this density does not span every irrep block, so the "
                                       "whole-system E cannot be built from it.");
    sys->AccumulateExchangeAll(X);
}

void Vxc::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
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
