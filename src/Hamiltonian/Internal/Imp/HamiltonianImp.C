// File: HamiltonianImp.C  General matrix implementation of a Hamiltonian operator.
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
module qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.stl_io;
import qchem.Streamable;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

template <class T> tHamiltonianImp<T>::tHamiltonianImp()
    : itsIsPolarized(false)
    , itsIsRelativistic(false)
{};

template <class T> void tHamiltonianImp<T>::Add(tStatic_HT<T>* p)
{
    itsSHTs.push_back(std::unique_ptr<tStatic_HT<T>>(p));
    itsIsPolarized    = itsIsPolarized    || p->IsPolarized();
    // AND, not OR: the virial holds only if EVERY term is Coulombic (V1.27).  One PP term kills it.
    itsIsVirialValid  = itsIsVirialValid  && p->IsVirialValid();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
    itsPreservesReal  = itsPreservesReal  && p->PreservesReal();   // AND: one SOC/A-field term flips all blocks complex
}
template <class T> void tHamiltonianImp<T>::Add(tDynamic_HT<T>* p)
{
    itsDHTs.push_back(std::unique_ptr<tDynamic_HT<T>>(p));
    itsIsPolarized    = itsIsPolarized    || p->IsPolarized();
    // AND, not OR: the virial holds only if EVERY term is Coulombic (V1.27).  One PP term kills it.
    itsIsVirialValid  = itsIsVirialValid  && p->IsVirialValid();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
    itsPreservesReal  = itsPreservesReal  && p->PreservesReal();   // AND: one SOC/A-field term flips all blocks complex
}
template <class T> void tHamiltonianImp<T>::Add(tDynamic_HF_HT<T>* p)
{
    itsHF_HTs.push_back(std::unique_ptr<tDynamic_HF_HT<T>>(p));
    itsIsPolarized    = itsIsPolarized    || p->IsPolarized();
    // AND, not OR: the virial holds only if EVERY term is Coulombic (V1.27).  One PP term kills it.
    itsIsVirialValid  = itsIsVirialValid  && p->IsVirialValid();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
    itsPreservesReal  = itsPreservesReal  && p->PreservesReal();   // AND: one SOC/A-field term flips all blocks complex
}

template <class T> hmat_t<T> tHamiltonianImp<T>::GetMatrix(const tobs_t<T>* bs,const Spin& S,const tChargeDensity<T>* cd,const tbs_t<T>* wholeBasis)
{
    // Layer-2 lineage guard: never build a Fock from a SUPERSEDED density (a previous iteration's, or a stale
    // copy).  The active (live-head) density is trivially active; a superseded one trips here at the exact
    // call site instead of silently returning a plausible-but-wrong matrix.  See ChargeDensity::Lineage.
    assert((!cd || cd->isActive()) && "Hamiltonian::GetMatrix computing with a superseded charge density");
    int n=bs->GetNumFunctions();
    hmat_t<T> H=blazem::zeroH<T>(n);
    for (auto& t:itsSHTs) H+=t->GetMatrix(bs,S);                       // static: no density
    // Leave these terms out if we don't have guess for the charge density.
    if (cd)
    {
        for (auto& t:itsDHTs)   H+=t->GetMatrix(bs,S,cd);             // per-irrep dynamic (DFT/fitted)
        for (auto& t:itsHF_HTs) H+=t->GetMatrix(bs,S,cd,wholeBasis);  // whole-system HF (needs the composite basis)
    }
    return H;
}


// The real-block fold (Step 3c-2): the same assembly as the native GetMatrix, through each term's
// Static/Dynamic_HT_RealBlock capability face (V1.6/V1.7 cross-cast; a term without the face is a wiring
// error and throws LOUDLY -- every periodic term has carried it since 3c-1).  Only ever called on the
// <dcmplx> instantiation (the real WF children of a complex run); molecular Hamiltonians never see it.
namespace
{
template <class HT> const Static_HT_RealBlock& StaticRealOf(const HT& t)
{
    auto* rb=dynamic_cast<const Static_HT_RealBlock*>(&t);
    if (!rb) throw std::logic_error("Hamiltonian real-block assembly: a static term does not carry the "
                                    "Static_HT_RealBlock capability (RealComplexPlan Step 3c-1)");
    return *rb;
}
template <class HT> const Dynamic_HT_RealBlock& DynamicRealOf(const HT& t)
{
    auto* rb=dynamic_cast<const Dynamic_HT_RealBlock*>(&t);
    if (!rb) throw std::logic_error("Hamiltonian real-block assembly: a dynamic term does not carry the "
                                    "Dynamic_HT_RealBlock capability (RealComplexPlan Step 3c-1)");
    return *rb;
}
} //anon

template <class T> hmat_t<double> tHamiltonianImp<T>::GetMatrix(const tobs_t<double>* bs,const Spin& S,
                                                                const tChargeDensity<dcmplx>* cd,
                                                                const tbs_t<dcmplx>*)
{
    assert((!cd || cd->isActive()) && "Hamiltonian::GetMatrix computing with a superseded charge density");
    int n=bs->GetNumFunctions();
    hmat_t<double> H=blazem::zeroH<double>(n);
    for (auto& t:itsSHTs) H+=StaticRealOf(*t).GetMatrix(bs,S);
    if (cd)
    {
        for (auto& t:itsDHTs) H+=DynamicRealOf(*t).GetMatrix(bs,S,cd);
        // No whole-system HF on the periodic path (exact exchange is real-molecular-only by construction);
        // a real block asked of an HF term would be a wiring error.
        if (!itsHF_HTs.empty())
            throw std::logic_error("Hamiltonian real-block assembly: whole-system HF terms cannot serve a "
                                   "real TRIM block (no periodic exact exchange)");
    }
    return H;
}

template <class T> EnergyBreakdown tHamiltonianImp<T>::GetTotalEnergy( const tDM_CD<T>* cd ) const
{
    assert(cd);
    assert(cd->isActive() && "Hamiltonian::GetTotalEnergy computing with a superseded charge density");
    EnergyBreakdown e;
    for (auto& t:itsSHTs)  t->GetEnergy(e,cd);
    for (auto& t:itsDHTs)  t->GetEnergy(e,cd);
    for (auto& t:itsHF_HTs)  t->GetEnergy(e,cd);
    return e;
}


template <class T> std::ostream& tHamiltonianImp<T>::Write(std::ostream& os) const
{
    if (itsIsPolarized) os << "Polarized ";
    if (itsIsRelativistic) os << "Relativistic ";
    os << "Hamiltonian with " << itsSHTs.size() << " static terms:" << std::endl;
    os << itsSHTs;
    os << "and " << itsDHTs.size() << " dynamic terms:" << std::endl;
    os << itsDHTs;
    os << "and " << itsHF_HTs.size() << " Hartree-Fock terms:" << std::endl;
    os << itsHF_HTs;
    return os;
}

template class tHamiltonianImp<double>;
template class tHamiltonianImp<dcmplx>;

} //namespace
