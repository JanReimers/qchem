// File: Hamiltonian/Internal/Imp/Vxc_Forces.C  The XC term's force at fixed D (doc/GeometryOptimizationPlan.md §13.2; lane B 2026-11-09).
// Its own TU so the force stays out of Imp/PWTerms_XC.C, which lane S's V-XC-ONE-TERM sprint edits.
//
// THE ENERGY.  E_xc = Σ_g w_g e(ρ↑(r_g), ρ↓(r_g)), the channel rasters the term's DensitySampler collocates from D, and the sampler's
// Matrix(bs, v) is the EXACT ADJOINT of that collocation (its class invariant: the two faces are one object).  So at fixed D only the
// basis functions move, and
//     dE_xc/dR|_D = Σ_g w_g Σ_σ v^σ(r_g) ∂ρ_σ(r_g)/∂R = Σ_blocks,σ Re Tr(D_σ ∂_R Matrix(bs, v^σ))   at the FIXED v^σ = ∂e/∂ρ_σ,
// the sampler's MakeMatrixPullback with the v the refresh built for D.  (On an atom-centred mesh the points and weights move too: that
// strategy refuses, so it is a MissingGradients() entry, never a silent zero.)
//
// THE ROUTES.  A plane-wave block on the ball-fit route: an EXACT ZERO (the fitter's overlap tensor has no centre to move).  The raw
// raster gather of a GPW block: refused (its jets are not built).  The δ quadrature (Becke or uniform mesh): refused.
module;
#include <stdexcept>
#include <string>
#include <type_traits>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.ChargeDensity;
import qchem.ChargeDensity.DensitySampler;
import qchem.Structure;

namespace qchem::Hamiltonian
{

namespace
{
//! The sampler's REAL-BLOCK face, for a real TRIM child of a complex run: the periodic sampler's (abstract -> abstract cross-cast).
const ChargeDensity::DensitySampler* PeriodicOf(const ChargeDensity::tDensitySampler<dcmplx>* s) {return dynamic_cast<const ChargeDensity::DensitySampler*>(s);}
//! The molecular path has no real TRIM child (its blocks are its own scalar): never asked, but the branch that asks must compile.
const ChargeDensity::DensitySampler* PeriodicOf(const ChargeDensity::tDensitySampler<double>*) {return nullptr;}
} // namespace

template <class T> std::string tVxc_Quadrature<T>::MissingGradient(const tbs_t<T>* bs, const Structure*) const
{
    if (!bs) return "Vxc_Quadrature: no basis to differentiate";
    std::string why;
    ForEachBlock(bs, [&](const auto* b)
    {
        if (!why.empty()) return;
        using BT=std::remove_const_t<std::remove_pointer_t<decltype(b)>>;
        bool can=false;
        if constexpr (std::is_same_v<BT,tobs_t<T>>) can=itsSampler->CanPullbackMatrix(b);
        else if (const auto* p=PeriodicOf(itsSampler.get())) can=p->CanPullbackMatrix(b);
        // The ROUTE refuses, not the block: the reason names no block, so a symmetry-blocked run lists the same entry as its unblocked twin.
        if (!can) why="Vxc_Quadrature: no geometry gradient (its XC quadrature's route gives no matrix pullback)";
    });
    return why;
}

template <class T> template <class U> rvec3vec_t tVxc_Quadrature<T>::PullbackBlock(const tobs_t<U>* bs, const Spin& s, const hmat_t<U>& D, const Structure* st) const
{
    const rvec_t& v=itsRefreshed.v[s==Spin::Down ? 1 : 0];   // a folded doublet (Spin::None) is the field either channel sees: Up (as MakeMatrix)
    if (v.size()==0) throw std::logic_error("tVxc_Quadrature::Pullback: the refresh built no v_xc for spin "+std::to_string(int(s))+" (AddForces checks the density first)");
    auto r=[&]
    {
        if constexpr (std::is_same_v<U,T>) return itsSampler->MakeMatrixPullback(bs, v, D, st);
        else
        {
            const auto* p=PeriodicOf(itsSampler.get());
            if (!p) throw std::logic_error("tVxc_Quadrature::Pullback: a real TRIM block needs the periodic sampler's real-block face");
            return p->MakeMatrixPullback(bs, v, D, st);
        }
    }();
    if (!r) throw std::logic_error("tVxc_Quadrature::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

template <class T> void tVxc_Quadrature<T>::AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const
{
    // A READ of what the refresh built: v^σ at THIS density (RefreshForDensity(bs, D) ran BuildBlocks), so the field is the derivative
    // of the energy at the D the force is taken at.
    if (!(DensityID(D)==itsRefreshed.id))
        throw std::logic_error("tVxc_Quadrature::AddForces: the refresh built v_xc for another density (refresh with the basis for D first)");
    const rvec3vec_t g=D->DM_Pullback(this, st);
    if (g.size()!=F.size()) throw std::logic_error("tVxc_Quadrature::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g[A];
}

// The members this TU defines (the class itself is instantiated in Imp/PWTerms_XC.C, where its other members live).
template std::string tVxc_Quadrature<double>::MissingGradient(const tbs_t<double>*, const Structure*) const;
template std::string tVxc_Quadrature<dcmplx>::MissingGradient(const tbs_t<dcmplx>*, const Structure*) const;
template void tVxc_Quadrature<double>::AddForces(rvec3vec_t&, const tDM_CD<double>*, const Structure*) const;
template void tVxc_Quadrature<dcmplx>::AddForces(rvec3vec_t&, const tDM_CD<dcmplx>*, const Structure*) const;
template rvec3vec_t tVxc_Quadrature<double>::PullbackBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;
template rvec3vec_t tVxc_Quadrature<dcmplx>::PullbackBlock<dcmplx>(const tobs_t<dcmplx>*, const Spin&, const hmat_t<dcmplx>&, const Structure*) const;
template rvec3vec_t tVxc_Quadrature<dcmplx>::PullbackBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;

} // namespace qchem::Hamiltonian
