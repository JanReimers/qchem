// File: HamiltonianImp.C  General matrix implementation of a Hamiltonian operator.
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <type_traits>   // is_same_v (MakeTermMatrices real-block twin)
module qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Hamiltonian.Internal.Term;     // DensityMatrixRequirement (folded at Add)
import qchem.Hamiltonian.Internal.ACBN0;    // ACBN0 + HubbardProjection (MakeACBN0Estimator)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.DensityMatrixDistance;   // DeltaRho without a quadrature
import qchem.stl_io;
import qchem.Streamable;
import qchem.Blaze;
import qchem.Parallel;   // OuterFor (InducedFock's coupled pairs, KAxisPlan P1)

namespace qchem::Hamiltonian
{

template <class T> tHamiltonianImp<T>::tHamiltonianImp(SpinGroup g)
    : itsSpinGroup(g)
    , itsIsRelativistic(false)
{};

template <class T> void tHamiltonianImp<T>::Add(tStatic_HT<T>* p)
{
    itsSHTs.push_back(std::unique_ptr<tStatic_HT<T>>(p));
    // AND, not OR: the virial holds only if EVERY term is Coulombic (V1.27).  One PP term kills it.
    itsIsVirialValid  = itsIsVirialValid  && p->IsVirialValid();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
    itsPreservesReal  = itsPreservesReal  && p->PreservesReal();   // AND: one SOC/A-field term flips all blocks complex
}
template <class T> std::unique_ptr<ACBN0Estimator> tHamiltonianImp<T>::MakeACBN0Estimator() const
{
    for (const auto& t : itsDHTs)
        if (auto* h=dynamic_cast<HubbardProjection*>(t.get())) return std::make_unique<ACBN0>(*h);
    return nullptr;
}

// Abstract -> abstract cross-cast: the +U term realises HubbardChannels; nothing here names Hubbard_U.
template <class T> const HubbardChannels* tHamiltonianImp<T>::GetHubbardChannels() const
{
    for (const auto& t : itsDHTs)
        if (auto* h=dynamic_cast<const HubbardChannels*>(t.get())) return h;
    return nullptr;
}
// The same cross-cast for the WRITE face (LinearResponsePlan §3d Q6).
template <class T> HubbardUTarget* tHamiltonianImp<T>::GetHubbardUTarget()
{
    for (const auto& t : itsDHTs)
        if (auto* h=dynamic_cast<HubbardUTarget*>(t.get())) return h;
    return nullptr;
}

namespace {

//! \brief The Hamiltonian's linearisation, folded over its dynamic terms' \c tResponse_HT capabilities
//! (doc/LinearResponsePlan.md H1/H2).  Non-owning: the terms belong to the Hamiltonian.
template <class T> class FoldedResponseKernel : public ResponseKernel<T>
{
public:
    FoldedResponseKernel(std::vector<const tResponse_HT<T>*> terms, const tbs_t<T>* wholeBasis,
                         const tChargeDensity<T>* D0)
        : itsTerms(std::move(terms)), itsWholeBasis(wholeBasis), itsD0(D0)
    {
        if (!itsWholeBasis) throw std::invalid_argument("ResponseKernel: the whole (composite) basis is required");
    }
    virtual std::unique_ptr<TransitionFock<T>> InducedFock(const TransitionDensity<T>& delta) const override
    {
        for (const auto* t : itsTerms) t->RefreshForDensity(itsWholeBasis, itsD0, delta);   // the response phase
        // The spin irreps δ resolves: {Up, Down}, or {None} -- the folded doublet (SpinIrrepsOf's rule).
        const std::vector<Spin> spins = (delta.Channel(Spin::Up) && delta.Channel(Spin::Down))
                                      ? std::vector<Spin>{Spin::Up, Spin::Down} : std::vector<Spin>{Spin::None};
        auto rule=delta.Coupling();
        auto dF=std::make_unique<AO_TransitionFock<T>>(rule);
        const auto blocks=itsWholeBasis->template Iterate<tobs_t<T>>();
        // THE COUPLED PAIRS RUN CONCURRENTLY (doc/KAxisPlan.md P1).  The refresh above filled every term's δ-derived state, so
        // GetMatrix(bra, ket, s, δ) is a READ (the concurrency contract on RefreshForDensity); each pair's δF goes to its own
        // slot and is added IN THE SERIAL ORDER below, so the result is bit-identical at any thread count.
        struct Pair { Spin s; const tobs_t<T>* bra; const tobs_t<T>* ket; };
        std::vector<Pair> pairs;
        for (const Spin& s : spins)
            for (const auto* ket : blocks)
                for (const auto* bra : blocks)
                    if (rule->Couples(bra->GetSymmetry(), ket->GetSymmetry())) pairs.push_back({s, bra, ket});
        std::vector<mat_t<T>> each(pairs.size());
        qchem::OuterFor(pairs.size(), [&](size_t i)
        {
            const Pair& p=pairs[i];
            mat_t<T> m(p.bra->GetNumFunctions(), p.ket->GetNumFunctions(), T(0));
            for (const auto* t : itsTerms) m+=t->GetMatrix(p.bra, p.ket, p.s, delta);
            each[i]=std::move(m);
        });
        for (size_t i=0;i<pairs.size();i++) dF->Add(pairs[i].bra->GetIrrep(pairs[i].s), pairs[i].ket->GetIrrep(pairs[i].s), each[i]);
        return dF;
    }
private:
    std::vector<const tResponse_HT<T>*> itsTerms;
    const tbs_t<T>*                     itsWholeBasis;
    const tChargeDensity<T>*            itsD0;
};

//! Collect \a t's response capability, or its description into \a missing.
template <class T, class Term> void CollectResponse(const Term& t, std::vector<const tResponse_HT<T>*>& out,
                                                   std::ostringstream& missing)
{
    if (auto* r=dynamic_cast<const tResponse_HT<T>*>(&t)) {out.push_back(r); return;}
    std::ostringstream name;
    t.Write(name);
    missing << "\n    " << name.str();
}

} // namespace

// Static terms are never asked: they do not depend on D.  Every DYNAMIC term must answer (ruling Q2).
template <class T> std::unique_ptr<ResponseKernel<T>> tHamiltonianImp<T>::MakeResponseKernel(const tbs_t<T>* wholeBasis,
                                                                                            const tChargeDensity<T>* D0) const
{
    std::vector<const tResponse_HT<T>*> terms;
    std::ostringstream missing;
    for (const auto& t : itsDHTs)   CollectResponse<T>(*t, terms, missing);
    for (const auto& t : itsHF_HTs) CollectResponse<T>(*t, terms, missing);
    if (!missing.str().empty())
        throw std::logic_error("Hamiltonian::MakeResponseKernel: these density-dependent terms cannot be linearised "
                               "yet (no tResponse_HT capability), so this Hamiltonian has no response kernel:"
                               + missing.str());
    return std::make_unique<FoldedResponseKernel<T>>(std::move(terms), wholeBasis, D0);
}

template <class T> void tHamiltonianImp<T>::Add(tDynamic_HT<T>* p)
{
    itsDHTs.push_back(std::unique_ptr<tDynamic_HT<T>>(p));
    // AND, not OR: the virial holds only if EVERY term is Coulombic (V1.27).  One PP term kills it.
    itsIsVirialValid  = itsIsVirialValid  && p->IsVirialValid();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
    itsPreservesReal  = itsPreservesReal  && p->PreservesReal();   // AND: one SOC/A-field term flips all blocks complex
    // OR: one term that needs D (the analytic (ij|δ) Hartree) makes the whole Hamiltonian need it.  Abstract -> abstract.
    if (auto* r=dynamic_cast<const DensityMatrixRequirement*>(p)) itsTermNeedsDensityMatrix = itsTermNeedsDensityMatrix || r->RequiresDensityMatrix();
}
template <class T> void tHamiltonianImp<T>::HoldNonlinearMesh(mesh_t m)
{
    if (!m) throw std::logic_error("tHamiltonianImp::HoldNonlinearMesh: no mesh (a factory that needs none holds none)");
    if (itsNonlinearMesh) throw std::logic_error("tHamiltonianImp::HoldNonlinearMesh: this Hamiltonian already holds a nonlinear integration "
                                                 "mesh -- one Hamiltonian, one mesh; its clients share it (V-NONLINEAR-MESH)");
    itsNonlinearMesh=std::move(m);
}

// THE OBSERVABLE (doc/OpenWork.md N1/T2; R1.0h): mu_A = Integral w_A (rho_up - rho_dn), in electrons, over the mesh's atom-centred
// partition.  FREE -- both channel rasters are cached for this density serial.  Empty when the run resolves no spin (m == 0
// identically), when the mesh has no site blocks, or when there is no mesh: the caller treats empty as "this run cannot answer",
// never as "zero".  It moved here from the XC term (V-NONLINEAR-MESH): the partition is the mesh's, not the functional's.
template <class T> rvec_t tHamiltonianImp<T>::SiteMoments(const tChargeDensity<T>* cd) const
{
    if (!itsNonlinearMesh || itsSpinGroup!=SpinGroup::Polarized) return rvec_t();
    const auto& m=*itsNonlinearMesh;
    return m.SiteIntegrals(rvec_t(m.RhoPol(cd,Spin::Up)-m.RhoPol(cd,Spin::Down)));
}

// D-MINDRHO (2026-10-22, user: "option A"): the L1 distance per electron on the run's nonlinear mesh -- the fraction of the electrons
// the step would still move, the same number for every mixer, representation and system (cached rasters: ~one extra raster per
// iteration), per spin channel on a polarized run (a moment still moving keeps the SCF unconverged).  With no mesh, the S-metric
// size of ΔD per electron over every block and spin, which needs no mesh at all.  The sampler caches ONE density serial (both
// channels of it), so BOTH of a's channels are copied before b's replace them -- two rasters per call, not four.
template <class T> double tHamiltonianImp<T>::DeltaRho(const tChargeDensity<T>* in, const tChargeDensity<T>* out) const
{
    double d=0.0;
    if (itsNonlinearMesh)
    {
        const auto& m=*itsNonlinearMesh;
        auto l1=[&](const rvec_t& fa, const rvec_t& fb)
        {
            rvec_t x(fa.size());
            for (size_t g=0; g<x.size(); g++) x[g]=std::fabs(fa[g]-fb[g]);
            return m.Integrate(x);
        };
        if (itsSpinGroup!=SpinGroup::Polarized)
        {
            const rvec_t fa=m.Rho(in);
            d=l1(fa, m.Rho(out));
        }
        else
        {
            const rvec_t au=m.RhoPol(in,Spin::Up), ad=m.RhoPol(in,Spin::Down);
            const rvec_t bu=m.RhoPol(out,Spin::Up), bd=m.RhoPol(out,Spin::Down);
            d=l1(au,bu)+l1(ad,bd);
        }
    }
    else
    {
        const auto* a=dynamic_cast<const ChargeDensity::tDensityMatrixDistance<T>*>(in);
        const auto* b=dynamic_cast<const ChargeDensity::tMixableDensity<T>*>(out);
        if (!a || !b) throw std::logic_error("tHamiltonianImp::DeltaRho: no nonlinear integration mesh in this Hamiltonian and no density "
                                             "matrix in these densities, so nothing can measure the change");
        d=std::sqrt(a->SquaredDistanceTo(*b));
    }
    const double N=out->GetTotalCharge();
    return N!=0.0 ? d/N : d;
}

template <class T> void tHamiltonianImp<T>::Add(tDynamic_HF_HT<T>* p)
{
    itsHF_HTs.push_back(std::unique_ptr<tDynamic_HF_HT<T>>(p));
    // AND, not OR: the virial holds only if EVERY term is Coulombic (V1.27).  One PP term kills it.
    itsIsVirialValid  = itsIsVirialValid  && p->IsVirialValid();
    itsIsRelativistic = itsIsRelativistic || p->IsRelativistic();
    itsPreservesReal  = itsPreservesReal  && p->PreservesReal();   // AND: one SOC/A-field term flips all blocks complex
    // A whole-system J/K term contracts D itself, so it needs D unless it says otherwise (DensityMatrixRequirement).  Abstract -> abstract.
    const auto* r=dynamic_cast<const DensityMatrixRequirement*>(p);
    itsTermNeedsDensityMatrix = itsTermNeedsDensityMatrix || (r ? r->RequiresDensityMatrix() : true);
}

template <class T> void tHamiltonianImp<T>::Add(tTerm<T>&& t)
{
    const int filled=int(bool(t.statik))+int(bool(t.dynamic))+int(bool(t.hf));
    if (filled!=1) throw std::logic_error("tHamiltonianImp::Add(Term): a term record holds "+std::to_string(filled)+" terms (exactly one expected)");
    if (t.mesh) HoldNonlinearMesh(std::move(t.mesh));   // the mesh BEFORE its first client term
    if (t.statik)  Add(t.statik.release());
    if (t.dynamic) Add(t.dynamic.release());
    if (t.hf)      Add(t.hf.release());
}

template <class T> hmat_t<T> tHamiltonianImp<T>::GetMatrix(const tobs_t<T>* bs,const Spin& S,const tChargeDensity<T>* cd,const tbs_t<T>* wholeBasis) const
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
        for (auto& t:itsDHTs)   H+=t->GetMatrix(bs,S,DensityID(cd));  // per-irrep dynamic (DFT/fitted)
        for (auto& t:itsHF_HTs) H+=t->GetMatrix(bs,S,DensityID(cd),wholeBasis);  // whole-system HF (needs the composite basis)
    }
    return H;
}


// The real-block fold (Step 3c-2): the same assembly as the native GetMatrix, through each term's
// real block face, tStatic_CC<double> / tDynamic_CC<double> (V1.6/V1.7 cross-cast; a term without the face is a wiring
// error and throws LOUDLY -- every periodic term has carried it since 3c-1).  Only ever called on the
// <dcmplx> instantiation (the real WF children of a complex run); molecular Hamiltonians never see it.
namespace
{
template <class HT> const tStatic_CC<double>& StaticRealOf(const HT& t)
{
    auto* rb=dynamic_cast<const tStatic_CC<double>*>(&t);   // abstract to abstract: the term's real block face (V1.35)
    if (!rb) throw std::logic_error("Hamiltonian real-block assembly: a static term does not carry the "
                                    "real block face tStatic_CC<double> (RealComplexPlan Step 3c-1)");
    return *rb;
}
template <class HT> const tDynamic_CC<double>& DynamicRealOf(const HT& t)
{
    auto* rb=dynamic_cast<const tDynamic_CC<double>*>(&t);   // abstract to abstract: the term's real block face (V1.35)
    if (!rb) throw std::logic_error("Hamiltonian real-block assembly: a dynamic term does not carry the "
                                    "real block face tDynamic_CC<double> (RealComplexPlan Step 3c-1)");
    return *rb;
}
} //anon

template <class Imp> hmat_t<double> HamiltonianImp_RealBlock<Imp>::GetMatrix(const tobs_t<double>* bs,const Spin& S,
                                                                             const tChargeDensity<dcmplx>* cd,
                                                                             const tbs_t<dcmplx>*) const
{
    const auto& itsSHTs=self().itsSHTs; const auto& itsDHTs=self().itsDHTs; const auto& itsHF_HTs=self().itsHF_HTs;
    assert((!cd || cd->isActive()) && "Hamiltonian::GetMatrix computing with a superseded charge density");
    int n=bs->GetNumFunctions();
    hmat_t<double> H=blazem::zeroH<double>(n);
    for (auto& t:itsSHTs) H+=StaticRealOf(*t).GetMatrix(bs,S);
    if (cd)
    {
        for (auto& t:itsDHTs) H+=DynamicRealOf(*t).GetMatrix(bs,S,DensityID(cd));
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
    RefreshForDensity(nullptr, cd);   // the contract (KAxisPlan P2): every term is refreshed for cd before it is asked
    EnergyBreakdown e;
    for (auto& t:itsSHTs)  t->GetEnergy(e,cd);
    for (auto& t:itsDHTs)  t->GetEnergy(e,cd);
    for (auto& t:itsHF_HTs)  t->GetEnergy(e,cd);
    // The integrated site moments are the MESH's observable (V-NONLINEAR-MESH), written with the energy so the trace and the
    // observer read them off the breakdown they already hold (R1.0h); free -- the energy pass cached both channel rasters.
    e.charge.siteMoments=SiteMoments(cd);
    return e;
}


// THE GRID AUDIT's per-term fold (doc/GridEpsilonPlan.md §3): GetTotalEnergy's / GetMatrix's loops, one entry per term instead of a sum.
namespace
{
template <class HT, class CD> TermEnergy EnergyOfTerm(const HT& t, const CD* cd)
{
    EnergyBreakdown e;
    t.GetEnergy(e,cd);
    TermEnergy r;
    for (const auto& [name,term] : e.Terms()) r.entries += (r.entries.empty() ? "" : "+") + name;
    r.E=e.GetTotalEnergy();
    std::ostringstream os; os<<t;                                          // the term's own description, first line, trimmed
    std::string l=os.str(); l=l.substr(0,l.find('\n'));
    const size_t a=l.find_first_not_of(" \t"); r.label = a==std::string::npos ? r.entries : l.substr(a);
    return r;
}
} //anon

template <class T> std::vector<TermEnergy> tHamiltonianImp<T>::MakeTermEnergies(const tDM_CD<T>* cd) const
{
    assert(cd && cd->isActive());
    // A whole-system HF term (HF, DHF, or a density model's exact J, coulomb.route Exact) has no raster in it, so a mesh audit has nothing
    // of it to compare: it is LEFT OUT, as the hub's MakeTermForces note says (its forces follow the energies' terms).  MakeTermMatrices
    // leaves it out too, so the two lists stay index for index.
    std::vector<TermEnergy> out;
    for (auto& t:itsSHTs) out.push_back(EnergyOfTerm(*t,cd));
    for (auto& t:itsDHTs) out.push_back(EnergyOfTerm(*t,cd));
    return out;
}

template <class T> std::vector<hmat_t<T>> tHamiltonianImp<T>::MakeTermMatrices(const tobs_t<T>* bs,const Spin& S,const tChargeDensity<T>* cd) const
{
    assert(cd && cd->isActive());
    std::vector<hmat_t<T>> out;   // the whole-system HF terms left out, as MakeTermEnergies does
    for (auto& t:itsSHTs) out.push_back(t->GetMatrix(bs,S));
    for (auto& t:itsDHTs) out.push_back(t->GetMatrix(bs,S,DensityID(cd)));
    return out;
}

template <class T> std::vector<hmat_t<double>> tHamiltonianImp<T>::MakeTermMatrices(const tobs_t<double>* bs,const Spin& S,
                                                                                   const tChargeDensity<dcmplx>* cd) const
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        assert(cd && cd->isActive());
        if (!itsHF_HTs.empty()) throw std::logic_error("Hamiltonian::MakeTermMatrices: whole-system HF terms cannot serve a real TRIM block");
        std::vector<hmat_t<double>> out;
        for (auto& t:itsSHTs) out.push_back(StaticRealOf(*t).GetMatrix(bs,S));
        for (auto& t:itsDHTs) out.push_back(DynamicRealOf(*t).GetMatrix(bs,S,DensityID(cd)));
        return out;
    }
    else throw std::logic_error("Hamiltonian::MakeTermMatrices: a real block of a REAL Hamiltonian is the native overload");
}

template <class T> std::ostream& tHamiltonianImp<T>::Write(std::ostream& os) const
{
    if (itsSpinGroup==SpinGroup::Polarized) os << "Polarized ";
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
template class HamiltonianImp_RealBlock<tHamiltonianImp<dcmplx>>;

} //namespace
