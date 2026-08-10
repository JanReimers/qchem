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

const rsmat_t& Dynamic_HF_HT_Imp::GetMatrix(const robs_t* bs,const Spin&,const rChargeDensity* cd,
                                            const rbs_t* wholeBasis) const
{
    // An HF term is inherently whole-system (canonical-pair ScatterBoth), so the composite basis is
    // required; there is no valid null-basis caller (the per-irrep fallback would fetch non-canonical blocks
    // and hit the §3c cache guard).  Fail loudly rather than route to a path that cannot work.
    if (!wholeBasis)
        throw std::runtime_error("HF term: the whole-system Fock build requires the composite basis "
                                 "(the cross-irrep view) -- GetMatrix was called with a null wholeBasis.");
    // Stash the run-stable whole basis on first use (GetEnergy has no basis argument and needs the same
    // banked contraction).  It must never CHANGE: itsJKs is keyed by BasisSetID and built by walking this
    // composite, so a second basis would have its blocks served from the first one's contraction -- a
    // plausible-looking Fock built from the wrong cross-irrep view.  THROW rather than assert, for the same
    // reason as R1.4: build/Release is -DNDEBUG, so an assert here would be a no-op in the configuration we
    // actually test, and the null check three lines up already throws on the same parameter.
    if (!itsWholeBasis) itsWholeBasis=wholeBasis;
    else if (itsWholeBasis!=wholeBasis)
        throw std::runtime_error("HF term: the whole-system (composite) basis changed mid-run.  This term "
                                 "latched a different basis on its first Fock build, and its per-irrep "
                                 "blocks (itsJKs) are keyed by BasisSetID against THAT basis -- serving them "
                                 "for a different composite would silently mix cross-irrep views.  A term "
                                 "belongs to one wavefunction; do not share it across two.");
    ContractAll(cd);
    return itsJKs.at(bs->BasisSetID());
}

void Dynamic_HF_HT_Imp::ContractAll(const rChargeDensity* cd) const
{
    assert(itsWholeBasis);
    if (cd->Version()==itsCD_Version && !itsJKs.empty()) return;   // already current for this density
    const rDM_CD* dm = dynamic_cast<const rDM_CD*>(cd);
    assert(dm && "HF term: density must be a rDM_CD");
    std::vector<const robs_t*> obs;                          // the whole basis's per-irrep blocks
    std::vector<rsmat_t>      X;                            // one zeroed block per irrep (same order as obs)
    for (auto* b:itsWholeBasis->Iterate<robs_t>())
    {
        obs.push_back(b);
        X.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    AccumulateAll(X,dm);                                    // canonical-pair scatter (diagonal + off-diagonal)
    const double w=Scale();
    itsJKs.clear();
    for (size_t k=0;k<obs.size();++k)
    {
        if (w!=1.0) X[k]*=w;                                // no-op multiply skipped for Coulomb (w==1)
        itsJKs[obs[k]->BasisSetID()]=std::move(X[k]);
    }
    itsCD_Version=cd->Version();
}

} //namespace
