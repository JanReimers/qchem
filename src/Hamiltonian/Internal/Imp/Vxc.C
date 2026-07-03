// File: Vxc.C  Hartree-Fock exchange potential
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  Let the charge density do the work.
//

// Whole-system RHF exchange (doc/ERI4Rework.md §5.4): the total density scatters itself across canonical
// irrep pairs (ScatterBoth on Exchange blocks), so K(j,i) is never built.  Blocks are stored already
// scaled by -1/2 (the RHF exchange coefficient), so GetMatrix can hand back a reference.
void Vxc::ContractAllExchange(const rChargeDensity* cd) const
{
    assert(itsWholeBasis);
    if (cd->Version()==itsAllVersion && !itsK.empty()) return;
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);
    assert(dm && "Vxc (HF exchange): density must be a DM_CD");
    std::vector<const obs_t*>   obs;
    std::vector<const ohfbs_t*> hf;
    std::vector<rsmat_t>        Kall;
    for (auto* b:itsWholeBasis->Iterate<obs_t>())
    {
        auto* h=dynamic_cast<const ohfbs_t*>(b);
        assert(h && "Vxc: irrep basis is not a Hartree-Fock orbital basis");
        obs.push_back(b);
        hf.push_back(h);
        Kall.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    dm->AccumulateExchangeAll(Kall,hf);
    itsK.clear();
    for (size_t k=0;k<obs.size();++k) { Kall[k]*=-0.5; itsK[obs[k]->BasisSetID()]=std::move(Kall[k]); }
    itsAllVersion=cd->Version();
}

const rsmat_t& Vxc::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd,const bs_t* wholeBasis) const
{
    // RHF exchange is whole-system (see Vee): the composite basis is required, no valid null-basis caller.
    if (!wholeBasis)
        throw std::runtime_error("Vxc (HF exchange): the whole-system Fock build requires the composite basis "
                                 "(the cross-irrep view) -- GetMatrix was called with a null wholeBasis.");
    if (!itsWholeBasis) itsWholeBasis=wholeBasis;
    ContractAllExchange(cd);
    return itsK.at(bs->BasisSetID());
}

void Vxc::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    // E_x = 1/2 Tr(D.K_scaled) from this term's own whole-system (already -1/2 scaled) exchange blocks.
    ContractAllExchange(cd);
    te.Exc+=0.5*cd->DM_ContractBlocks(itsK);
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    os << "    Hartee-Fock exchange potential phi(r_1)*phi(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
