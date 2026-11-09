// File: Hamiltonian/Internal/Imp/Ven_Forces.C  The nuclear-attraction term's force (doc/GeometryOptimizationPlan.md §13.2; lane B
// 2026-11-08).  Its own TU so the force work stays out of Imp/Ven.C and the term files lane S's sprint edits.
//
// E_en = Tr(D V) at fixed D, and V moves with BOTH the basis centres and the nuclei, so ∂E_en/∂R_A = Tr(D ∂V/∂R_A) with both motions: each
// irrep block's from its basis face (Integrals_Nuclear::MakeNuclearPullback), summed by the density walk (tDM_CD::DM_Pullback).
module;
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Structure;

namespace qchem::Hamiltonian
{

std::string Ven::MissingGradient(const tbs_t<double>* bs, const Structure*) const
{
    if (!bs) return "Ven: no basis to differentiate";
    for (const auto* b : bs->template Iterate<robs_t>())
        if (!b->CanPullbackNuclear()) return "Ven: the basis block "+b->Name()+" gives no nuclear pullback";
    return "";
}

rvec3vec_t Ven::Pullback(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure* st) const
{
    const auto r=bs->MakeNuclearPullback(D, st);
    if (!r) throw std::logic_error("Ven::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

void Ven::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const
{
    const rvec3vec_t g=D->DM_Pullback(this, st);
    if (g.size()!=F.size()) throw std::logic_error("Ven::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g[A];
}

} // namespace qchem::Hamiltonian
