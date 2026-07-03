// File: Vee.C  Electron-Electron Coulomb potential
module;
#include <cassert>
#include <iostream>
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
//  This is where we calculate the overlap of the fit basis functions with
//  the real coulomb potential:
//              /
//  Vreal(r1) = | dr2 Ro(r2)/r12 .
//              /
//  Where ro is the charge density.
//

// Assemble the WHOLE-system Coulomb once per density (doc/ERI4Rework.md §5.4).  The context supplies every
// ab-irrep basis; the density scatters itself across canonical irrep pairs (ScatterBoth), so only the
// canonical ERI4 blocks are ever built/cached -- halving Coulomb ERI RAM+build vs the old per-irrep pass,
// which fetched both J(i,j) and J(j,i).  Result is bit-close (~1e-13) to the old Fock, not bit-identical:
// the old J(j,i) was an independent build, here it is J(i,j)^T.
void Vee::ContractAllDirect(const rChargeDensity* cd) const
{
    assert(itsWholeBasis);
    if (cd->Version()==itsAllVersion && !itsJ.empty()) return;   // already current for this density
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);
    assert(dm && "Vee (HF Coulomb): density must be a DM_CD");
    std::vector<const obs_t*>   obs;                         // the whole basis's per-irrep blocks
    std::vector<const ohfbs_t*> hf;                          // HF face, for the ERI scatter
    std::vector<rsmat_t>        Jall;                        // one zeroed Coulomb block per irrep
    for (auto* b:itsWholeBasis->Iterate<obs_t>())
    {
        auto* h=dynamic_cast<const ohfbs_t*>(b);
        assert(h && "Vee: irrep basis is not a Hartree-Fock orbital basis");
        obs.push_back(b);
        hf.push_back(h);
        Jall.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    dm->AccumulateDirectAll(Jall,hf);                        // canonical-pair ScatterBoth into every block
    itsJ.clear();
    for (size_t k=0;k<obs.size();++k) itsJ[obs[k]->BasisSetID()]=std::move(Jall[k]);
    itsAllVersion=cd->Version();
}

const rsmat_t& Vee::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd,const bs_t* wholeBasis) const
{
    // HF Coulomb is inherently whole-system (canonical-pair ScatterBoth), so the composite basis is
    // required; there is no valid null-basis caller (the per-irrep fallback would fetch non-canonical
    // blocks and hit the §3c cache guard).  Fail loudly rather than route to a path that cannot work.
    if (!wholeBasis)
        throw std::runtime_error("Vee (HF Coulomb): the whole-system Fock build requires the composite basis "
                                 "(the cross-irrep view) -- GetMatrix was called with a null wholeBasis.");
    if (!itsWholeBasis) itsWholeBasis=wholeBasis;                 // stash the run-stable whole basis
    ContractAllDirect(cd);
    return itsJ.at(bs->BasisSetID());
}

void Vee::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    // E_ee = 1/2 Tr(D.J), taken from THIS term's own whole-system Coulomb blocks -- no per-irrep GetMatrix
    // round-trip through DM_Contract (which is what kept Vee tied to the 3-arg GetMatrix / tDynamic_CC).
    ContractAllDirect(cd);
    te.Eee=0.5*cd->DM_ContractBlocks(itsJ);
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

std::ostream& Vee::Write(std::ostream& os) const
{
    os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
