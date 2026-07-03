// File: HF_HT.C  Shared whole-system machinery for the 4-index Hartree-Fock terms (Vee, Vxc).
//
// The version guard + composite-basis walk + per-irrep block cache used to be copy-pasted between
// Vee::ContractAllDirect and Vxc::ContractAllExchange (they differed by exactly one contraction call and
// Vxc's scale).  It now lives once here on Dynamic_HF_HT_Imp; each term supplies only AccumulateAll (the
// Direct-vs-Exchange line) and, optionally, Scale.  See doc/ERI4Rework.md §5.4.
module;
#include <cassert>
#include <cstddef>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

const rsmat_t& Dynamic_HF_HT_Imp::GetMatrix(const obs_t* bs,const Spin&,const rChargeDensity* cd,
                                            const bs_t* wholeBasis) const
{
    // An HF term is inherently whole-system (canonical-pair ScatterBoth), so the composite basis is
    // required; there is no valid null-basis caller (the per-irrep fallback would fetch non-canonical blocks
    // and hit the §3c cache guard).  Fail loudly rather than route to a path that cannot work.
    if (!wholeBasis)
        throw std::runtime_error("HF term: the whole-system Fock build requires the composite basis "
                                 "(the cross-irrep view) -- GetMatrix was called with a null wholeBasis.");
    if (!itsWholeBasis) itsWholeBasis=wholeBasis;                 // stash the run-stable whole basis
    ContractAll(cd);
    return itsBlocks.at(bs->BasisSetID());
}

void Dynamic_HF_HT_Imp::ContractAll(const rChargeDensity* cd) const
{
    assert(itsWholeBasis);
    if (cd->Version()==itsAllVersion && !itsBlocks.empty()) return;   // already current for this density
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);
    assert(dm && "HF term: density must be a DM_CD");
    std::vector<const obs_t*>   obs;                         // the whole basis's per-irrep blocks
    std::vector<const ohfbs_t*> hf;                          // HF face, for the ERI scatter
    std::vector<rsmat_t>        X;                           // one zeroed block per irrep
    for (auto* b:itsWholeBasis->Iterate<obs_t>())
    {
        auto* h=dynamic_cast<const ohfbs_t*>(b);
        assert(h && "HF term: irrep basis is not a Hartree-Fock orbital basis");
        obs.push_back(b);
        hf.push_back(h);
        X.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    AccumulateAll(X,hf,dm);                                  // canonical-pair ScatterBoth into every block
    const double w=Scale();
    itsBlocks.clear();
    for (size_t k=0;k<obs.size();++k)
    {
        if (w!=1.0) X[k]*=w;                                // no-op multiply skipped for Coulomb (w==1)
        itsBlocks[obs[k]->BasisSetID()]=std::move(X[k]);
    }
    itsAllVersion=cd->Version();
}

} //namespace
