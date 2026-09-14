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

const rsmat_t& Dynamic_HF_HT_Imp::GetMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd,
                                            const rbs_t* wholeBasis) const
{
    if (!wholeBasis)
        throw std::runtime_error("HF term: the whole-system Fock build requires the composite basis "
                                 "(the cross-irrep view) -- GetMatrix was called with a null wholeBasis.");
    if (!itsWholeBasis) itsWholeBasis=wholeBasis;
    else if (itsWholeBasis!=wholeBasis)
        throw std::runtime_error("HF term: the whole-system (composite) basis changed mid-run.  This term "
                                 "latched a different basis on its first Fock build, and its per-irrep "
                                 "blocks (itsJKs) are keyed by BasisSetID against THAT basis -- serving them "
                                 "for a different composite would silently mix cross-irrep views.  A term "
                                 "belongs to one wavefunction; do not share it across two.");
    return ContractAll(cd, CacheSpin(s)).at(bs->BasisSetID());
}

// One (density, spin) contraction: the canonical-pair scatter over the whole basis, scaled, keyed per irrep.
// \a s is a CacheSpin: None for a spin-blind operator (one set of blocks for the run), Up/Down for a
// per-channel one -- so a polarized exchange keeps K_up and K_dn side by side under one density serial.
const std::map<std::string,rsmat_t>& Dynamic_HF_HT_Imp::ContractAll(const rChargeDensity* cd, const Spin& s) const
{
    assert(itsWholeBasis);
    Blocks& b=itsJKs[s];
    if (cd->Version()==b.version && !b.jk.empty()) return b.jk;   // already current for this density
    const rDM_CD* dm=DensityFor(cd, s);
    std::vector<const robs_t*> obs;                          // the whole basis's per-irrep blocks
    std::vector<rsmat_t>      X;                            // one zeroed block per irrep (same order as obs)
    for (auto* ob:itsWholeBasis->Iterate<robs_t>())
    {
        obs.push_back(ob);
        X.push_back(blazem::zero<double>(ob->GetNumFunctions()));
    }
    AccumulateAll(X,dm);                                    // canonical-pair scatter (diagonal + off-diagonal)
    const double w=Scale(s);
    b.jk.clear();
    for (size_t k=0;k<obs.size();++k)
    {
        if (w!=1.0) X[k]*=w;                                // no-op multiply skipped for Coulomb (w==1)
        b.jk[obs[k]->BasisSetID()]=std::move(X[k]);
    }
    b.version=cd->Version();
    return b.jk;
}

} //namespace
