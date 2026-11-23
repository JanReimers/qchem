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
// raster gather of a GPW block: the gather's own pullback at the fixed field (PairProjection::applyRawAdjointPullback, lane B 2026-11-19).
// The δ quadrature on a FREE atom-centred (Becke) mesh: per block the functions moving and the points riding their owners, plus the mesh's
// weight response contracted against e (lane B 2026-11-19, READY 4).  An imposed-symmetry mesh refuses by name.
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

template <class T> std::string tVxc_DFT<T>::MissingGradient(const tbs_t<T>* bs, const Structure*) const
{
    if (!bs) return "cVxc_DFT: no basis to differentiate";
    std::string why;
    ForEachBlock(bs, [&](const auto* b)
    {
        if (!why.empty()) return;
        using BT=std::remove_const_t<std::remove_pointer_t<decltype(b)>>;
        bool can=false;
        if constexpr (std::is_same_v<BT,tobs_t<T>>) can=itsSampler->CanPullbackMatrix(b);
        else if (const auto* p=PeriodicOf(itsSampler.get())) can=p->CanPullbackMatrix(b);
        // The ROUTE refuses, not the block: the reason names no block, so a symmetry-blocked run lists the same entry as its unblocked twin.
        if (!can) why="cVxc_DFT: no geometry gradient (its XC quadrature gives no matrix pullback on this run: an IMPOSED-symmetry Becke mesh, whose orbit-rotated grids and star-averaged rho have no derivative yet, or a route without one)";
    });
    return why;
}

template <class T> template <class U> rvec3vec_t tVxc_DFT<T>::PullbackBlock(const tobs_t<U>* bs, const Spin& s, const hmat_t<U>& D, const Structure* st) const
{
    const rvec_t& v=itsRefreshed.v[s==Spin::Down ? 1 : 0];   // a folded doublet (Spin::None) is the field either channel sees: Up (as MakeMatrix)
    if (v.size()==0) throw std::logic_error("tVxc_DFT::Pullback: the refresh built no v_xc for spin "+std::to_string(int(s))+" (AddForces checks the density first)");
    auto r=[&]
    {
        if constexpr (std::is_same_v<U,T>) return itsSampler->MakeMatrixPullback(bs, v, D, st);
        else
        {
            const auto* p=PeriodicOf(itsSampler.get());
            if (!p) throw std::logic_error("tVxc_DFT::Pullback: a real TRIM block needs the periodic sampler's real-block face");
            return p->MakeMatrixPullback(bs, v, D, st);
        }
    }();
    if (!r) throw std::logic_error("tVxc_DFT::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

template <class T> const rvec_t& tVxc_DFT<T>::FieldFor(const Spin& s, const char* who) const
{
    const rvec_t& v=itsRefreshed.v[s==Spin::Down ? 1 : 0];   // a folded doublet (Spin::None) is the field either channel sees: Up (as MakeMatrix)
    if (v.size()==0) throw std::logic_error(std::string("tVxc_DFT::")+who+": the refresh built no v_xc for spin "+std::to_string(int(s)));
    return v;
}
template <class T> ChargeDensity::SampledBlocks<T> tVxc_DFT<T>::Sampled(const ChargeDensity::DensityBlocks<T>& b, const char* who) const
{
    ChargeDensity::SampledBlocks<T> s;
    for (const auto& x : b.blocks) s.blocks.push_back({x.bs, &FieldFor(x.s, who), x.D});
    for (const auto& x : b.real)   s.real  .push_back({x.bs, &FieldFor(x.s, who), x.D});
    return s;
}
template <class T> rvec3vec_t tVxc_DFT<T>::PullbackBlocks(const ChargeDensity::DensityBlocks<T>& b, const Structure* st) const
{
    const auto r=itsSampler->MakeMatrixPullbackBlocks(Sampled(b, "PullbackBlocks"), st);
    if (!r) throw std::logic_error("tVxc_DFT::PullbackBlocks: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

template <class T> void tVxc_DFT<T>::AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const
{
    // A READ of what the refresh built: v^σ at THIS density (RefreshForDensity(bs, D) ran BuildBlocks), so the field is the derivative
    // of the energy at the D the force is taken at.
    if (!(DensityID(D)==itsRefreshed.id))
        throw std::logic_error("tVxc_DFT::AddForces: the refresh built v_xc for another density (refresh with the basis for D first)");
    const rvec3vec_t g=D->DM_Pullback(this, st);
    if (g.size()!=F.size()) throw std::logic_error("tVxc_DFT::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g[A];
    // An ATOM-CENTRED mesh's weights move with the atoms: their response, contracted against e (the blocks above held the weights; lane B
    // 2026-11-19, the controller-approved split).  A mesh whose weights cannot respond refuses its blocks' pullbacks, so it never reaches here.
    if (itsSampler->HasQuadratureResponse())
    {
        const auto q=itsSampler->MakeQuadraturePullback(EnergyDensityAt(D), st);
        if (!q) throw std::logic_error("tVxc_DFT::AddForces: "+q.Error().why+" (MissingGradients() names it before a force is asked)");
        if (q->gradient.size()!=F.size()) throw std::logic_error("tVxc_DFT::AddForces: the weight response and the force array differ in atom count");
        for (size_t A=0;A<F.size();A++) F[A]-=q->gradient[A];
    }
}

// The members this TU defines (the class itself is instantiated in Imp/PWTerms_XC.C, where its other members live).
template std::string tVxc_DFT<double>::MissingGradient(const tbs_t<double>*, const Structure*) const;
template std::string tVxc_DFT<dcmplx>::MissingGradient(const tbs_t<dcmplx>*, const Structure*) const;
template void tVxc_DFT<double>::AddForces(rvec3vec_t&, const tDM_CD<double>*, const Structure*) const;
template void tVxc_DFT<dcmplx>::AddForces(rvec3vec_t&, const tDM_CD<dcmplx>*, const Structure*) const;
template const rvec_t& tVxc_DFT<double>::FieldFor(const Spin&, const char*) const;
template const rvec_t& tVxc_DFT<dcmplx>::FieldFor(const Spin&, const char*) const;
template ChargeDensity::SampledBlocks<double> tVxc_DFT<double>::Sampled(const ChargeDensity::DensityBlocks<double>&, const char*) const;
template ChargeDensity::SampledBlocks<dcmplx> tVxc_DFT<dcmplx>::Sampled(const ChargeDensity::DensityBlocks<dcmplx>&, const char*) const;
template rvec3vec_t tVxc_DFT<double>::PullbackBlocks(const ChargeDensity::DensityBlocks<double>&, const Structure*) const;
template rvec3vec_t tVxc_DFT<dcmplx>::PullbackBlocks(const ChargeDensity::DensityBlocks<dcmplx>&, const Structure*) const;
template rvec3vec_t tVxc_DFT<double>::PullbackBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;
template rvec3vec_t tVxc_DFT<dcmplx>::PullbackBlock<dcmplx>(const tobs_t<dcmplx>*, const Spin&, const hmat_t<dcmplx>&, const Structure*) const;
template rvec3vec_t tVxc_DFT<dcmplx>::PullbackBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;

} // namespace qchem::Hamiltonian
