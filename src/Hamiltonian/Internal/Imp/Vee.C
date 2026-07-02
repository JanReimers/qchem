// File: Vee.C  Electron-Electron Coulomb potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <string>

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

rsmat_t Vee::CalcMatrix(const obs_t* bs,const Spin&,const rChargeDensity* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    auto hf_bs = dynamic_cast<const ohfbs_t*>(bs);
    assert(hf_bs);
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);   // HF J needs the density matrix (not a fit seed)
    assert(dm && "Vee (HF Coulomb): density must be a DM_CD");
    rsmat_t Jab=blazem::zero<double>(bs->GetNumFunctions());
    dm->AccumulateDirect(Jab,hf_bs);
    return Jab;
}

// Assemble the WHOLE-system Coulomb once per density (doc/ERI4Rework.md §5.4).  The context supplies every
// ab-irrep basis; the density scatters itself across canonical irrep pairs (ScatterBoth), so only the
// canonical ERI4 blocks are ever built/cached -- halving Coulomb ERI RAM+build vs the old per-irrep pass,
// which fetched both J(i,j) and J(j,i).  Result is bit-close (~1e-13) to the old Fock, not bit-identical:
// the old J(j,i) was an independent build, here it is J(i,j)^T.
void Vee::EnsureWholeSystem(const rChargeDensity* cd) const
{
    assert(!itsBases.empty());
    if (cd->Version()==itsAllVersion && !itsJ.empty()) return;   // already current for this density
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);
    assert(dm && "Vee (HF Coulomb): density must be a DM_CD");
    const size_t N=itsBases.size();
    std::vector<const ohfbs_t*> hf; hf.reserve(N);           // HF face, for the ERI scatter
    std::vector<rsmat_t>        Jall; Jall.reserve(N);       // one zeroed Coulomb block per irrep
    for (auto* b:itsBases)
    {
        auto* h=dynamic_cast<const ohfbs_t*>(b);
        assert(h && "Vee: context irrep basis is not a Hartree-Fock orbital basis");
        hf.push_back(h);
        Jall.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    dm->AccumulateDirectAll(Jall,hf);                        // canonical-pair ScatterBoth into every block
    itsJ.clear();
    for (size_t k=0;k<N;++k) itsJ[itsBases[k]->BasisSetID()]=std::move(Jall[k]);
    itsAllVersion=cd->Version();
}

const rsmat_t& Vee::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd,const HamiltonianContext& ctx) const
{
    if (ctx.irrepBases.empty()) return Dynamic_HT_Imp::GetMatrix(bs,s,cd);   // no cross-irrep view: per-irrep path
    if (itsBases.empty()) itsBases=ctx.irrepBases;                           // stash the run-stable ab-irrep list
    newCD(cd);
    EnsureWholeSystem(cd);
    return itsJ.at(bs->BasisSetID());
}

const rsmat_t& Vee::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    // Energy / other context-free callers.  Once a Fock build has stashed the irrep list, this density
    // (typically the post-diagonalization density the energy is evaluated on) gets the SAME symmetry-banked
    // whole-system build -- so the non-canonical ERI4 blocks are never materialized here either.  Before
    // any context has arrived (stand-alone tests), fall back to the per-irrep base path.
    if (itsBases.empty()) return Dynamic_HT_Imp::GetMatrix(bs,s,cd);
    EnsureWholeSystem(cd);
    return itsJ.at(bs->BasisSetID());
}

void Vee::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    te.Eee=0.5*cd->DM_Contract(this,cd);
    te.EeeFit    = 0.0;
    te.EeeFitFit = 0.0;
}

std::ostream& Vee::Write(std::ostream& os) const
{
    os << "    Coulomb potential ro(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
