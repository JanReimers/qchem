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

// Same-spin: the block's OWN channel, or the whole (folded) density for a Spin::None block.  A polarized
// block handed a density that does not resolve its spin is a composition error, not a case.
const rDM_CD* Vxc::DensityFor(const rChargeDensity* cd, const Spin& s) const
{
    const rDM_CD* dm=ChargeDensity::DM_ChannelOf(cd,s);
    if (!dm) throw std::runtime_error("HF exchange: a spin-polarized block asked for its channel of a density "
                                      "that does not resolve spin (or carries no density matrix).");
    return dm;
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    os << "    Hartree-Fock exchange potential phi(r_1)*phi(r_2)/r_12 (same-spin: -K[D_sigma])" << std::endl;
    return os;
}

} //namespace
