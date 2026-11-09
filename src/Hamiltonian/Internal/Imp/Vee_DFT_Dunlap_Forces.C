// File: Hamiltonian/Internal/Imp/Vee_DFT_Dunlap_Forces.C  The Dunlap Hartree term's force (doc/GeometryOptimizationPlan.md §13.2; lane B
// 2026-11-08).  Its own TU so the force work stays out of Imp/Vee_DFT.C and the term files lane S's sprint edits.
//
// THE ENERGY.  E = Tr(D V) − ½ cᵀGc = c·b − ½ cᵀGc, with b_k = Σ_ij D_ij (ij|k), G = (k|l) the Coulomb metric, and c the charge-constrained
// fit (Dunlap-Connolly-Sabin): stationary in c under g·c = N, g_k = ∫f_k.
//
// THE FORCE.  dE/dR|_D = ∂E/∂R|_c + (∂E/∂c)·dc/dR.  Stationarity makes ∂E/∂c = λ g with g fixed (a function's integral does not move with
// it), and g·c = N, so the second term is λ dN/dR: the fit's response enters only through its charge.  At fixed D that charge,
// N = Tr(DS), moves with the basis: λ Tr(D ∂S/∂R).  It is LEFT OUT, by the contract on tGeometryGradient_HT::AddForces: this term's
// matrix is V = Σ c_k (ij|k), not its D-gradient V + λS (λS shifts every level by λ and moves no orbital), so the Pulay sum's W lacks
// the λ D that would cancel it.  What remains is
//     Σ_ij D_ij c_k ∂(ij|k)/∂R_A − ½ Σ_kl c_k c_l ∂(k|l)/∂R_A,
// the first the density walk of the fitter's repulsion pullback (orbital AND fit centres moving), the second its self-repulsion pullback.
// So the per-term fixed-D finite difference is this plus λ Tr(D ∂S/∂R), and the gate is the whole force against the re-converged energy
// (Calculation/tests/ForceFacade.C, Xα at α = 0: Dunlap is the only electron-electron term).
module;
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.FittedCD;
import qchem.Structure;

namespace qchem::Hamiltonian
{

std::string Vee_DFT_Dunlap::MissingGradient(const tbs_t<double>* bs, const Structure*) const
{
    if (!bs) return "Vee_DFT_Dunlap: no basis to differentiate";
    for (const auto* b : bs->template Iterate<robs_t>())
    {
        const auto* dft=dynamic_cast<const odftbs_t*>(b);
        if (!dft) return "Vee_DFT_Dunlap: the basis block "+b->Name()+" has no DFT tier";
        if (!itsFittedChargeDensity->CanPullbackRepulsion(dft))
            return "Vee_DFT_Dunlap: the basis block "+b->Name()+" gives no 3-centre repulsion pullback against its fit";
    }
    if (!itsFittedChargeDensity->CanPullbackSelfRepulsion()) return "Vee_DFT_Dunlap: the fit basis gives no Coulomb-metric pullback";
    return "";
}

rvec3vec_t Vee_DFT_Dunlap::Pullback(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure* st) const
{
    const auto* dft=dynamic_cast<const odftbs_t*>(bs);
    if (!dft) throw std::logic_error("Vee_DFT_Dunlap::Pullback: the block has no DFT tier (MissingGradients() names it before a force is asked)");
    const auto r=itsFittedChargeDensity->MakeRepulsionPullback(dft, D, st);
    if (!r) throw std::logic_error("Vee_DFT_Dunlap::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

void Vee_DFT_Dunlap::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const
{
    RequireFitFor(DensityID(D));   // a READ: c is the fit the refresh made of this D (tHamiltonianAt ran it)
    const rvec3vec_t g3=D->DM_Pullback(this, st);                        // Σ D c ∂(ij|k)/∂R
    const auto g2=itsFittedChargeDensity->MakeSelfRepulsionPullback(st);  // ½ Σ c c ∂(k|l)/∂R
    if (!g2) throw std::logic_error("Vee_DFT_Dunlap::AddForces: "+g2.Error().why+" (MissingGradients() names it before a force is asked)");
    if (g3.size()!=F.size() || g2->size()!=F.size())
        throw std::logic_error("Vee_DFT_Dunlap::AddForces: the pullbacks and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g3[A]-(*g2)[A];
}

} // namespace qchem::Hamiltonian
