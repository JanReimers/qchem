// File: Vxc.C  Hartree-Fock exchange potential
module;
#include <iostream>
#include <cassert>
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
//  Let the charge density do the work.
//

rsmat_t Vxc::CalcMatrix(const obs_t* bs,const Spin&,const rChargeDensity* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    auto hf_bs = dynamic_cast<const ohfbs_t*>(bs);
    assert(hf_bs);
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);   // HF K needs the density matrix (not a fit seed)
    assert(dm && "Vxc (HF exchange): density must be a DM_CD");
    rsmat_t Kab=blazem::zero<double>(bs->GetNumFunctions());
    dm->AccumulateExchange(Kab,hf_bs);
    return Kab*-0.5;
}

// Whole-system RHF exchange (doc/ERI4Rework.md §5.4): the total density scatters itself across canonical
// irrep pairs (ScatterBoth on Exchange blocks), so K(j,i) is never built.  Blocks are stored already
// scaled by -1/2 (the RHF exchange coefficient), so GetMatrix can hand back a reference.
void Vxc::EnsureWholeSystem(const rChargeDensity* cd) const
{
    assert(!itsBases.empty());
    if (cd->Version()==itsAllVersion && !itsK.empty()) return;
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);
    assert(dm && "Vxc (HF exchange): density must be a DM_CD");
    const size_t N=itsBases.size();
    std::vector<const ohfbs_t*> hf; hf.reserve(N);
    std::vector<rsmat_t>        Kall; Kall.reserve(N);
    for (auto* b:itsBases)
    {
        auto* h=dynamic_cast<const ohfbs_t*>(b);
        assert(h && "Vxc: context irrep basis is not a Hartree-Fock orbital basis");
        hf.push_back(h);
        Kall.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    dm->AccumulateExchangeAll(Kall,hf);
    itsK.clear();
    for (size_t k=0;k<N;++k) { Kall[k]*=-0.5; itsK[itsBases[k]->BasisSetID()]=std::move(Kall[k]); }
    itsAllVersion=cd->Version();
}

const rsmat_t& Vxc::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd,const HamiltonianContext& ctx) const
{
    if (ctx.irrepBases.empty()) return Dynamic_HT_Imp::GetMatrix(bs,s,cd);   // no cross-irrep view
    if (itsBases.empty()) itsBases=ctx.irrepBases;
    newCD(cd);
    EnsureWholeSystem(cd);
    return itsK.at(bs->BasisSetID());
}

const rsmat_t& Vxc::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    if (itsBases.empty()) return Dynamic_HT_Imp::GetMatrix(bs,s,cd);         // energy path before any Fock build
    EnsureWholeSystem(cd);
    return itsK.at(bs->BasisSetID());
}
void Vxc::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    te.Exc+=0.5*cd->DM_Contract(this,cd);
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    os << "    Hartee-Fock exchange potential phi(r_1)*phi(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
