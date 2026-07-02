// File: VxcPol.C  Polarized HF exchange potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.Symmetry.Spin;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real exchange potential,  Vxc(ro(r)), where ro is the charge density.
//
// The Hamiltonain matrix elements are calculated
//             /
//  Vxc(i,j) = | dr Vxcfit(ro(r)) Oi(r) Oj(r) .
//             /
//
//           = Sum  { Ck <Oi|Vk|Oj> } .
//
//  This last part is carried out by the base class FitImplementation.
rsmat_t VxcPol::CalcMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    if  (s==Spin::None)
    {
        std::cerr << "PolarizedHartreeFockVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    auto hf_bs = dynamic_cast<const ohfbs_t*>(bs);
    assert(hf_bs);

    const Polarized_CD* PolExactCD =  dynamic_cast<const Polarized_CD*>(cd);
    assert(PolExactCD);
    const DM_CD* SpinCD   = PolExactCD->GetChargeDensity(s); //Get CD for this spin direction
    rsmat_t Kab=blazem::zero<double>(bs->GetNumFunctions());
    SpinCD->AccumulateExchange(Kab,hf_bs);
    return Kab*-1.0;
}

// Whole-system UHF exchange for ONE spin: that spin's density scatters itself across canonical irrep pairs
// (ScatterBoth on Exchange blocks), so K(j,i) is never built.  Cached per (spin,irrep), already scaled by
// -1.  The version clear drops BOTH spins when the density changes; each spin builds on first request.
void VxcPol::EnsureWholeSystem(const rChargeDensity* cd, const Spin& s) const
{
    assert(!itsBases.empty());
    if (cd->Version()!=itsAllVersion) { itsK.clear(); itsAllVersion=cd->Version(); }
    if (itsK.count(s)) return;                                   // this spin already current
    const Polarized_CD* pcd = dynamic_cast<const Polarized_CD*>(cd);
    assert(pcd && "VxcPol: density must be polarized");
    const DM_CD* SpinCD = pcd->GetChargeDensity(s);
    const size_t N=itsBases.size();
    std::vector<const ohfbs_t*> hf; hf.reserve(N);
    std::vector<rsmat_t>        Kall; Kall.reserve(N);
    for (auto* b:itsBases)
    {
        auto* h=dynamic_cast<const ohfbs_t*>(b);
        assert(h && "VxcPol: context irrep basis is not a Hartree-Fock orbital basis");
        hf.push_back(h);
        Kall.push_back(blazem::zero<double>(b->GetNumFunctions()));
    }
    SpinCD->AccumulateExchangeAll(Kall,hf);
    auto& m=itsK[s];
    for (size_t k=0;k<N;++k) { Kall[k]*=-1.0; m[itsBases[k]->BasisSetID()]=std::move(Kall[k]); }
}

const rsmat_t& VxcPol::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd,const HamiltonianContext& ctx) const
{
    if (s==Spin::None) { std::cerr << "VxcPol::GetMatrix: unpolarized spin in a polarized term" << std::endl; exit(-1); }
    if (ctx.irrepBases.empty()) return Dynamic_HT_Imp_NoCache::GetMatrix(bs,s,cd);
    if (itsBases.empty()) itsBases=ctx.irrepBases;
    EnsureWholeSystem(cd,s);
    return itsK.at(s).at(bs->BasisSetID());
}

const rsmat_t& VxcPol::GetMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    if (itsBases.empty() || s==Spin::None) return Dynamic_HT_Imp_NoCache::GetMatrix(bs,s,cd);
    EnsureWholeSystem(cd,s);
    return itsK.at(s).at(bs->BasisSetID());
}
void VxcPol::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Exc += 0.5*cd->DM_Contract(this,cd); //This should sum K^alpha and K^beta.
}

std::ostream& VxcPol::Write(std::ostream& os) const
{
    os << "    Polarized Hartee-Fock exchange potential." << std::endl;
    return os;
}

} //namespace
