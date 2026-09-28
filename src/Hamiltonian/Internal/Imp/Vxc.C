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
void Vxc::AccumulateAll(std::vector<rsmat_t>& X,const ChargeDensity::tHF_System_CD<double>& sweep) const
{
    sweep.AccumulateExchangeAll(X);
}

// E_x = 1/2 Sum_sigma Tr(D_sigma.K_sigma_scaled) over the spin irreps the density resolves: {Up, Down} of a
// polarized composite (each channel against its own -K[D_sigma]), or {None} -- the folded doublet against
// -1/2 K[D_tot], which is the RHF energy.  One expression, no second type (V1.37 step 3).
void Vxc::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    double trDK=0.0;
    for (const Spin& s : ChargeDensity::SpinIrrepsOf(cd))
        trDK+=DensityFor(cd,s)->DM_ContractBlocks(ContractAll(cd,s));
    te.Add("Exc", 0.5*trDK, EnergyRole::Potential, trDK);   // quadratic: E = 1/2 Tr(D K), Tr(D V) = Tr(D K)
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    os << "    Hartree-Fock exchange potential phi(r_1)*phi(r_2)/r_12 (same-spin: -K[D_sigma])" << std::endl;
    return os;
}

} //namespace
