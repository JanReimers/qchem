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

// One (density, spin) contraction: the density's CacheSpin channel, presented through its HF sweep face.
// \a s is a CacheSpin: None for a spin-blind operator (one set of blocks for the run), Up/Down for a
// per-channel one -- so a polarized exchange keeps K_up and K_dn side by side under one density serial.
const std::map<std::string,rsmat_t>& Dynamic_HF_HT_Imp::ContractAll(const rChargeDensity* cd, const Spin& s) const
{
    Blocks& b=itsJKs[s];
    if (cd->Version()==b.version && !b.jk.empty()) return b.jk;   // already current: skip the channel lookup too
    // V1.6: ask for the WHOLE-SYSTEM exact-exchange face rather than calling through the general density
    // face.  Before, a density that could not span the irreps (a bare leaf) hit an assert-only default --
    // a silent NO-OP under -DNDEBUG, i.e. a zeroed J/K and a wrong Fock in the build we ship.
    auto* sys=dynamic_cast<const ChargeDensity::tHF_System_CD<double>*>(DensityFor(cd, s));
    if (!sys) throw std::runtime_error("HF term: this density does not span every irrep block, so the "
                                       "whole-system D cannot be built from it.");
    return Contract(b, *sys, cd->Version(), s);
}

const std::map<std::string,rsmat_t>& Dynamic_HF_HT_Imp::Contract(Blocks& b, const ChargeDensity::tHF_System_CD<double>& sweep,
                                                                 size_t version, const Spin& s) const
{
    assert(itsWholeBasis);
    if (version==b.version && !b.jk.empty()) return b.jk;   // already current for this operand
    std::vector<const robs_t*> obs;                          // the whole basis's per-irrep blocks
    std::vector<rsmat_t>      X;                            // one zeroed block per irrep (same order as obs)
    for (auto* ob:itsWholeBasis->Iterate<robs_t>())
    {
        obs.push_back(ob);
        X.push_back(blazem::zero<double>(ob->GetNumFunctions()));
    }
    AccumulateAll(X,sweep);                                 // canonical-pair scatter (diagonal + off-diagonal)
    const double w=Scale(s);
    b.jk.clear();
    for (size_t k=0;k<obs.size();++k)
    {
        if (w!=1.0) X[k]*=w;                                // no-op multiply skipped for Coulomb (w==1)
        b.jk[obs[k]->BasisSetID()]=std::move(X[k]);
    }
    b.version=version;
    return b.jk;
}

// The spin rule of BOTH terms, from their one answer (CacheSpin): Coulomb (None) reads the whole density,
// exchange (s) its own channel -- ChannelOf answers the whole density for None, so one line covers both.
// A polarized block handed a density that does not resolve its spin is a composition error, not a case.
const rDM_CD* Dynamic_HF_HT_Imp::DensityFor(const rChargeDensity* cd, const Spin& s) const
{
    const rDM_CD* dm=ChargeDensity::DM_ChannelOf(cd, CacheSpin(s));
    if (!dm) throw std::runtime_error("HF term: a spin-polarized block asked for its channel of a density that "
                                      "does not resolve spin (or the density carries no density matrix).");
    return dm;
}

} //namespace
