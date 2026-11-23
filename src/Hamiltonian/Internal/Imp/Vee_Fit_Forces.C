// File: Hamiltonian/Internal/Imp/Vee_Fit_Forces.C  The first-order fitted Hartree term's force (doc/GeometryOptimizationPlan.md §13.2;
// lane B 2026-11-09).  Its own TU so the force stays out of Imp/Vee_Fit.C, which lane S's sprint edits.
//
// THE ENERGY.  On an orthonormal fit basis (plane waves: the GPW and PW Hartree) the fit IS the projection, ρ̃ = P D, and
// E_H = ½ ρ̃ᵀ J ρ̃ with J the basis's diagonal Coulomb metric.  E_H is quadratic in D, and its matrix is the physics round
// V = P^† J ρ̃ (Repulsion).  At fixed D only P moves with the atoms (a plane wave has no centre, so J does not):
//     dE/dR|_D = ρ̃ᵀ J ∂P D = Re Tr(D ∂_R V[ρ̃ fixed]),
// the pullback of the physics round at the fixed fit, summed over the blocks, ONCE (the ½ and the two factors of P cancel).  So the
// term's force is -Σ_blocks MakeRepulsionPullback, and its per-term fixed-D finite difference is exactly that.
//
// THE ROUTES.  GPW: the Coulomb tensor's applyAdjointPullback, the jet walk of the gather (LatticeCollocation::IntegratePotentialPullback).
// A plane-wave block: exactly zero (no centre).  A δ-mesh fit (the molecular first-order Hartree): refused, its mesh rides with the atoms
// (V-NONLINEAR-MESH), so it is a MissingGradients() entry, never a silent zero.
module;
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Structure;
import qchem.BasisSet.Orbital_DFT_IBS;   // the block's DFT face (what the fitter's pullback takes)
import qchem.BasisSet.PairProjector;     // FoldBlockRef / FoldBlockD (the whole-density pullback's list)

namespace qchem::Hamiltonian
{

namespace
{
//! The block's DFT face on this run's fit axis \a T, for either block scalar \a U (abstract -> abstract).
template <class T, class U> const BasisSet::Orbital_DFT_IBS<U,T>* DFTof(const tobs_t<U>* b) {return dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,T>*>(b);}
} // namespace

template <class T> std::string Vee_Fit<T>::MissingGradient(const tbs_t<T>* bs, const Structure*) const
{
    if (!bs) return "Vee_Fit: no basis to differentiate";
    std::string why;
    ForEachBlock(bs, [&](const auto* b)
    {
        if (!why.empty()) return;
        const auto* dft=DFTof<T>(b);
        if (!dft) why="Vee_Fit: the basis block "+b->Name()+" has no DFT tier";
        else if (!itsFit->CanPullbackRepulsion(*dft)) why="Vee_Fit: the basis block "+b->Name()+" gives no repulsion pullback against its fit";
    });
    return why;
}

template <class T> template <class U> rvec3vec_t Vee_Fit<T>::PullbackBlock(const tobs_t<U>* bs, const Spin&, const hmat_t<U>& D, const Structure* st) const
{
    const auto* dft=DFTof<T>(bs);
    if (!dft) throw std::logic_error("Vee_Fit::Pullback: the block has no DFT tier (MissingGradients() names it before a force is asked)");
    const auto r=itsFit->MakeRepulsionPullback(*dft, D, st);
    if (!r) throw std::logic_error("Vee_Fit::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

// The fitter's view of a whole density: each block's DFT face and its D, on the face its scalar has.
namespace
{
template <class T> void FitBlocks(const ChargeDensity::DensityBlocks<T>& b, std::vector<BasisSet::FoldBlockRef>& refs,
                                  std::vector<BasisSet::FoldBlockD>& Ds, const char* who)
{
    for (const auto& x : b.blocks)
    {
        const auto* dft=DFTof<T>(x.bs);
        if (!dft) throw std::logic_error(std::string(who)+": the block has no DFT tier (MissingGradients() names it before a force is asked)");
        BasisSet::FoldBlockRef r; BasisSet::FoldBlockD d;
        if constexpr (std::is_same_v<T,dcmplx>) { r.bloch=dft; d.bloch=x.D; } else { r.real=nullptr; }
        refs.push_back(r); Ds.push_back(d);
    }
    for (const auto& x : b.real)
    {
        const auto* dft=DFTof<T>(x.bs);
        if (!dft) throw std::logic_error(std::string(who)+": the block has no DFT tier (MissingGradients() names it before a force is asked)");
        BasisSet::FoldBlockRef r; BasisSet::FoldBlockD d;
        r.real=dft; d.real=x.D;
        refs.push_back(r); Ds.push_back(d);
    }
}
} // namespace

template <class T> rvec3vec_t Vee_Fit<T>::PullbackBlocks(const ChargeDensity::DensityBlocks<T>& b, const Structure* st) const
{
    if constexpr (!std::is_same_v<T,dcmplx>) return tStatic_PC<T>::PullbackBlocks(b, st);
    else
    {
        std::vector<BasisSet::FoldBlockRef> refs; std::vector<BasisSet::FoldBlockD> Ds;
        FitBlocks(b, refs, Ds, "Vee_Fit::PullbackBlocks");
        const auto r=itsFit->MakeRepulsionPullbackBlocks(refs, Ds, st);
        if (!r) throw std::logic_error("Vee_Fit::PullbackBlocks: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
        return *r;
    }
}
template <class T> rmat3d_t Vee_Fit<T>::VirialBlocks(const ChargeDensity::DensityBlocks<T>& b, const Structure* st) const
{
    if constexpr (!std::is_same_v<T,dcmplx>) return tStatic_VC<T>::VirialBlocks(b, st);
    else
    {
        std::vector<BasisSet::FoldBlockRef> refs; std::vector<BasisSet::FoldBlockD> Ds;
        FitBlocks(b, refs, Ds, "Vee_Fit::VirialBlocks");
        const auto r=itsFit->MakeRepulsionVirialBlocks(refs, Ds, st);
        if (!r) throw std::logic_error("Vee_Fit::VirialBlocks: "+r.Error().why+" (MissingVirials() names it before a stress is asked)");
        return *r;
    }
}

template <class T> void Vee_Fit<T>::AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const
{
    this->RequireFitFor(DensityID(D));   // a READ: ρ̃ is the fit the refresh made of this D (tHamiltonianAt ran it)
    const rvec3vec_t g=D->DM_Pullback(this, st);
    if (g.size()!=F.size()) throw std::logic_error("Vee_Fit::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g[A];
}

// The members this TU defines (the class itself is instantiated in Imp/Vee_Fit.C, where its other members live).
template std::string Vee_Fit<double>::MissingGradient(const tbs_t<double>*, const Structure*) const;
template std::string Vee_Fit<dcmplx>::MissingGradient(const tbs_t<dcmplx>*, const Structure*) const;
template void Vee_Fit<double>::AddForces(rvec3vec_t&, const tDM_CD<double>*, const Structure*) const;
template void Vee_Fit<dcmplx>::AddForces(rvec3vec_t&, const tDM_CD<dcmplx>*, const Structure*) const;
template rvec3vec_t Vee_Fit<double>::PullbackBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;
template rvec3vec_t Vee_Fit<dcmplx>::PullbackBlock<dcmplx>(const tobs_t<dcmplx>*, const Spin&, const hmat_t<dcmplx>&, const Structure*) const;
template rvec3vec_t Vee_Fit<dcmplx>::PullbackBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;
template rvec3vec_t Vee_Fit<double>::PullbackBlocks(const ChargeDensity::DensityBlocks<double>&, const Structure*) const;
template rvec3vec_t Vee_Fit<dcmplx>::PullbackBlocks(const ChargeDensity::DensityBlocks<dcmplx>&, const Structure*) const;
template rmat3d_t   Vee_Fit<double>::VirialBlocks  (const ChargeDensity::DensityBlocks<double>&, const Structure*) const;
template rmat3d_t   Vee_Fit<dcmplx>::VirialBlocks  (const ChargeDensity::DensityBlocks<dcmplx>&, const Structure*) const;

} // namespace qchem::Hamiltonian
