// File: Hamiltonian/Internal/Imp/HamiltonianForces.C  The Hamiltonian's force, folded over its terms (doc/GeometryOptimizationPlan.md §13,
// §13.2; lane B 2026-11-08).  Its own TU so the force work stays out of the term files lane S's sprint edits.
//
// THE SUM.  F_A = −dE/dR_A at the converged density = Σ_terms −∂E_term/∂R_A |_D (each term's tGeometryGradient_HT, reached by the
// sanctioned abstract -> abstract cross-cast, as MakeResponseKernel reaches tResponse_HT) + Tr(W ∂S/∂R_A) (the orthonormality, or Pulay,
// term: the orbitals stay S-orthonormal as the basis moves; W = Σ f ε c cᵀ).  MissingGradients() names every term without the face, every
// face that refuses, and a basis without the overlap pullback the Pulay sum needs: never a silent zero.
module;
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Hamiltonian.Internal.GeometryGradient;
import qchem.ChargeDensity;
import qchem.Structure;

namespace qchem::Hamiltonian
{
namespace
{
//! The Pulay sum's per-block client: the overlap's pullback, asked of the block the density walk hands it.
template <class T> class OverlapPullback : public virtual tStatic_PC<T>
{
public:
    rvec3vec_t Pullback(const tobs_t<T>* bs, const Spin&, const hmat_t<T>& W, const Structure* st) const override
    {
        const auto r=bs->MakeOverlapPullback(W, st);
        if (!r) throw std::logic_error("Hamiltonian::Forces, the Pulay term: "+r.Error().why);
        return *r;
    }
};
//! A term's name: the first non-blank line of its Write, trimmed.
template <class Term> std::string NameOf(const Term& t)
{
    std::ostringstream os; t.Write(os);
    std::istringstream in(os.str());
    for (std::string line; std::getline(in, line); )
        if (const auto b=line.find_first_not_of(" \t"); b!=std::string::npos) return line.substr(b);
    return "(unnamed term)";
}
} // namespace

template <class T> std::vector<std::string> tHamiltonianImp<T>::MissingGradients(const tbs_t<T>* bs, const Structure* st) const
{
    std::vector<std::string> out;
    auto ask=[&](const auto& t)
    {
        if (const auto* g=dynamic_cast<const tGeometryGradient_HT<T>*>(&t))
        {
            if (std::string m=g->MissingGradient(bs, st); !m.empty()) out.push_back(m);
        }
        else out.push_back(NameOf(t)+": no geometry gradient (the term has no tGeometryGradient_HT face yet)");
    };
    for (const auto& t : itsSHTs)   ask(*t);
    for (const auto& t : itsDHTs)   ask(*t);
    for (const auto& t : itsHF_HTs) ask(*t);
    if (!bs) out.push_back("Pulay, Tr(W dS/dR): no basis to differentiate");
    else for (const auto* b : bs->template Iterate<tobs_t<T>>())
        if (!b->CanPullbackOverlap()) out.push_back("Pulay, Tr(W dS/dR): the basis block "+b->Name()+" gives no overlap pullback");
    return out;
}

template <class T> rvec3vec_t tHamiltonianImp<T>::MakeForces(const tbs_t<T>* bs, const tDM_CD<T>* D, const tDM_CD<T>* W,
                                                             const Structure* st) const
{
    if (!st || !D || !W) throw std::invalid_argument("Hamiltonian::Forces: needs the structure, the density matrix and W");
    if (const auto missing=MissingGradients(bs, st); !missing.empty())
    {
        std::string why;
        for (const auto& m : missing) why+="\n    "+m;
        throw std::logic_error("Hamiltonian::Forces: asked of a Hamiltonian that cannot give one (MissingGradients() says so):"+why);
    }
    rvec3vec_t F(st->GetNumAtoms(), rvec3_t(0,0,0));
    auto add=[&](const auto& t){ dynamic_cast<const tGeometryGradient_HT<T>&>(t).AddForces(F, D, st); };
    for (const auto& t : itsSHTs)   add(*t);
    for (const auto& t : itsDHTs)   add(*t);
    for (const auto& t : itsHF_HTs) add(*t);
    const OverlapPullback<T> pulay;
    const rvec3vec_t g=W->DM_Pullback(&pulay, st);
    if (g.size()!=F.size()) throw std::logic_error("Hamiltonian::Forces: the Pulay sum and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]+=g[A];
    return F;
}

template <class T> std::vector<TermForce> tHamiltonianImp<T>::MakeTermForces(const tbs_t<T>* bs, const tDM_CD<T>* D, const Structure* st) const
{
    if (!st || !D) throw std::invalid_argument("Hamiltonian::MakeTermForces: needs the structure and the density matrix");
    std::vector<TermForce> out;
    auto one=[&](const auto& t)
    {
        TermForce f{NameOf(t), "", {}};
        const auto* g=dynamic_cast<const tGeometryGradient_HT<T>*>(&t);
        f.missing = g ? g->MissingGradient(bs, st) : f.label+": no geometry gradient (the term has no tGeometryGradient_HT face yet)";
        if (f.missing.empty()) { f.F=rvec3vec_t(st->GetNumAtoms(), rvec3_t(0,0,0)); g->AddForces(f.F, D, st); }
        out.push_back(std::move(f));
    };
    for (const auto& t : itsSHTs)   one(*t);   // MakeTermEnergies' order: static, then dynamic ...
    for (const auto& t : itsDHTs)   one(*t);
    for (const auto& t : itsHF_HTs) one(*t);   // ... then the whole-system HF terms, which MakeTermEnergies leaves out (its grid audit has no raster in them)
    return out;
}

// Only THIS file's members: the class itself is explicitly instantiated in Imp/HamiltonianImp.C (a second class-level instantiation
// would be a duplicate definition).
template std::vector<std::string> tHamiltonianImp<double>::MissingGradients(const tbs_t<double>*, const Structure*) const;
template std::vector<std::string> tHamiltonianImp<dcmplx>::MissingGradients(const tbs_t<dcmplx>*, const Structure*) const;
template rvec3vec_t tHamiltonianImp<double>::MakeForces(const tbs_t<double>*, const tDM_CD<double>*, const tDM_CD<double>*, const Structure*) const;
template rvec3vec_t tHamiltonianImp<dcmplx>::MakeForces(const tbs_t<dcmplx>*, const tDM_CD<dcmplx>*, const tDM_CD<dcmplx>*, const Structure*) const;
template std::vector<TermForce> tHamiltonianImp<double>::MakeTermForces(const tbs_t<double>*, const tDM_CD<double>*, const Structure*) const;
template std::vector<TermForce> tHamiltonianImp<dcmplx>::MakeTermForces(const tbs_t<dcmplx>*, const tDM_CD<dcmplx>*, const Structure*) const;

} // namespace qchem::Hamiltonian
