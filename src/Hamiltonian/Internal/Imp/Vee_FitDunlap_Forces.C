// File: Hamiltonian/Internal/Imp/Vee_FitDunlap_Forces.C  The Dunlap Hartree term's force (doc/GeometryOptimizationPlan.md §13.2; lane B
// 2026-11-08; carried onto the V-ONE-OF-N 3b/3c shape by lane S: the term holds a Fitting::DunlapDensityFit, whose fitter owns c).
//
// THE ENERGY.  E = Tr(D V) − ½ cᵀGc = c·b − ½ cᵀGc, with b_k = Σ_ij D_ij (ij|k), G = (k|l) the Coulomb metric, and c the charge-constrained
// fit (Dunlap-Connolly-Sabin): stationary in c under g·c = N, g_k = ∫f_k, with N the FIXED electron count (user ruling 2026-11-08).
//
// THE FORCE.  dE/dR|_D = ∂E/∂R|_c + (∂E/∂c)·dc/dR.  Stationarity makes ∂E/∂c = λ g with g fixed (a function's integral does not move with
// it), and g·c = N, so the second term is λ dN/dR = 0: N is fixed, not Tr(DS).  What remains is
//     Σ_ij D_ij c_k ∂(ij|k)/∂R_A − ½ Σ_kl c_k c_l ∂(k|l)/∂R_A,
// the first the density walk of the fitter's repulsion pullback (orbital AND fit centres moving), the second its self-repulsion pullback.
// (Before the ruling the constraint was Tr(DS), and λ Tr(D ∂S/∂R) was the contract's stated exception; lane B retires that text.)
// So the per-term fixed-D finite difference is exactly this (lane B re-enables that gate), and the whole-force gate is against the re-converged energy
// (Calculation/tests/ForceFacade.C, Xα at α = 0: Dunlap is the only electron-electron term).
module;
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Structure;
import qchem.BasisSet.Orbital_DFT_IBS;   // the block's DFT face (what the fitter's pullback takes)

namespace qchem::Hamiltonian
{

std::string Vee_FitDunlap::MissingGradient(const tbs_t<double>* bs, const Structure*) const
{
    if (!bs) return "Vee_FitDunlap: no basis to differentiate";
    for (const auto* b : bs->template Iterate<robs_t>())
    {
        const auto* dft=dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,double>*>(b);
        if (!dft) return "Vee_FitDunlap: the basis block "+b->Name()+" has no DFT tier";
        if (!itsFit->CanPullbackRepulsion(*dft))
            return "Vee_FitDunlap: the basis block "+b->Name()+" gives no 3-centre repulsion pullback against its fit";
    }
    if (!itsFit->CanPullbackSelfRepulsion()) return "Vee_FitDunlap: the fit basis gives no Coulomb-metric pullback";
    return "";
}

rvec3vec_t Vee_FitDunlap::Pullback(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure* st) const
{
    const auto* dft=dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,double>*>(bs);
    if (!dft) throw std::logic_error("Vee_FitDunlap::Pullback: the block has no DFT tier (MissingGradients() names it before a force is asked)");
    const auto r=itsFit->MakeRepulsionPullback(*dft, D, st);
    if (!r) throw std::logic_error("Vee_FitDunlap::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

void Vee_FitDunlap::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const
{
    RequireFitFor(DensityID(D));   // a READ: c is the fit the refresh made of this D (tHamiltonianAt ran it)
    const rvec3vec_t g3=D->DM_Pullback(this, st);                        // Σ D c ∂(ij|k)/∂R
    const auto g2=itsFit->MakeSelfRepulsionPullback(st);                  // ½ Σ c c ∂(k|l)/∂R (the fitter's, scaled as SelfRepulsion)
    if (!g2) throw std::logic_error("Vee_FitDunlap::AddForces: "+g2.Error().why+" (MissingGradients() names it before a force is asked)");
    if (g3.size()!=F.size() || g2->size()!=F.size())
        throw std::logic_error("Vee_FitDunlap::AddForces: the pullbacks and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g3[A]-(*g2)[A];
}

} // namespace qchem::Hamiltonian
