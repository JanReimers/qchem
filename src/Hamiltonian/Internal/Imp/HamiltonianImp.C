// File: HamiltonianImp.C  General matrix implementation of a Hamiltonian operator.
module;
#include <cassert>
#include <iostream>
#include <memory>
module qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.IonIon;   // IonIon<double> (the ion-ion energy term)
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
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
}
template <class T> void tHamiltonianImp<T>::Add(tDynamic_HT<T>* p)
{
    itsDHTs.push_back(std::unique_ptr<tDynamic_HT<T>>(p));
    itsIsPolarized    = itsIsPolarized    || p->IsPolarized();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
}
template <class T> void tHamiltonianImp<T>::Add(tDynamic_HF_HT<T>* p)
{
    itsHF_HTs.push_back(std::unique_ptr<tDynamic_HF_HT<T>>(p));
    itsIsPolarized    = itsIsPolarized    || p->IsPolarized();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
}

// The molecular standard terms (Kinetic/Vnn/Ven) are double-only; the complex (plane-wave) Hamiltonian
// builds its terms explicitly, so this is NA there.
template <> void tHamiltonianImp<double>::InsertStandardTerms(const st_t & st)
{
    Add(new Kinetic);
    Add(new IonIon<double>(st));
    Add(new Ven(st));
}
template <> void tHamiltonianImp<dcmplx>::InsertStandardTerms(const st_t &)
{
    assert(false && "InsertStandardTerms: the complex Hamiltonian assembles its terms explicitly");
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
