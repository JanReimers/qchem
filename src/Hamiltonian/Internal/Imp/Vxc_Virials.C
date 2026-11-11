// File: Hamiltonian/Internal/Imp/Vxc_Virials.C  The XC term's VIRIAL at fixed D (doc/Records/CellHalf_Design.md §5; user ruling §6.2,
// 2026-11-10; lane H).  Its own TU, beside Imp/Vxc_Forces.C, so the strain work stays out of Imp/PWTerms_XC.C (lane S's sprint).
//
// THE ENERGY.  E_xc = Σ_g w_g e(ρ↑_g, ρ↓_g) on the sampler's channel rasters.  Under a homogeneous strain at fixed D:
//   - the QUADRATURE side (the pair counts held): on a raster fixed in fractional coordinates with weights Ω/N, every ρ_g scales as 1/Ω,
//     so ∂E/∂ε_ab = δ_ab (E_xc − Σ_σ Σ_g w_g v^σ_g ρ^σ_g), v^σ = ∂e/∂ρ_σ -- the pure dilation (a sampler whose points move refuses);
//   - the ORBITAL side: Σ_g w_g v_g ∂ρ_g at the fixed raster and v = Σ_blocks,σ the sampler's MakeMatrixVirial (its Matrix is the exact
//     adjoint of its Rho) -- a plane-wave block's is an exact zero, its pair counts being fixed by its G indices.
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
const ChargeDensity::DensitySampler* PeriodicOf(const ChargeDensity::tDensitySampler<dcmplx>* s) {return dynamic_cast<const ChargeDensity::DensitySampler*>(s);}
const ChargeDensity::DensitySampler* PeriodicOf(const ChargeDensity::tDensitySampler<double>*) {return nullptr;}
} // namespace

template <class T> std::string tVxc_Quadrature<T>::MissingVirial(const tbs_t<T>* bs, const Structure*) const
{
    if (!bs) return "Vxc_Quadrature: no basis to differentiate";
    if (!itsSampler->HasDilationVirial()) return "Vxc_Quadrature: no virial (its XC quadrature's points do not stay fractional under a strain)";
    std::string why;
    ForEachBlock(bs, [&](const auto* b)
    {
        if (!why.empty()) return;
        using BT=std::remove_const_t<std::remove_pointer_t<decltype(b)>>;
        bool can=false;
        if constexpr (std::is_same_v<BT,tobs_t<T>>) can=itsSampler->CanVirialMatrix(b);
        else if (const auto* p=PeriodicOf(itsSampler.get())) can=p->CanVirialMatrix(b);
        if (!can) why="Vxc_Quadrature: no virial (its XC quadrature's route gives no matrix virial)";
    });
    return why;
}

template <class T> template <class U> rmat3d_t tVxc_Quadrature<T>::VirialBlock(const tobs_t<U>* bs, const Spin& s, const hmat_t<U>& D, const Structure* st) const
{
    const rvec_t& v=itsRefreshed.v[s==Spin::Down ? 1 : 0];   // a folded doublet (Spin::None) is the field either channel sees: Up (as MakeMatrix)
    if (v.size()==0) throw std::logic_error("tVxc_Quadrature::Virial: the refresh built no v_xc for spin "+std::to_string(int(s))+" (AddVirial checks the density first)");
    auto r=[&]
    {
        if constexpr (std::is_same_v<U,T>) return itsSampler->MakeMatrixVirial(bs, v, D, st);
        else
        {
            const auto* p=PeriodicOf(itsSampler.get());
            if (!p) throw std::logic_error("tVxc_Quadrature::Virial: a real TRIM block needs the periodic sampler's real-block face");
            return p->MakeMatrixVirial(bs, v, D, st);
        }
    }();
    if (!r) throw std::logic_error("tVxc_Quadrature::Virial: "+r.Error().why+" (MissingVirials() names it before a stress is asked)");
    return *r;
}

template <class T> void tVxc_Quadrature<T>::AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const
{
    if (!(DensityID(D)==itsRefreshed.id))
        throw std::logic_error("tVxc_Quadrature::AddVirial: the refresh built v_xc for another density (refresh with the basis for D first)");
    if (!itsSampler->HasDilationVirial())
        throw std::logic_error("tVxc_Quadrature::AddVirial: the quadrature's points do not stay fractional (MissingVirials() names it first)");
    // The quadrature side: E_xc and Σ_σ ∫ v^σ ρ^σ on the SAME channel rasters GetEnergy integrates.
    rvec_t half;
    const Rasters r=ChannelRasters(D, half);
    const rvec_t& vUp=itsRefreshed.v[0];
    const rvec_t& vDn=itsRefreshed.v[1].size() ? itsRefreshed.v[1] : itsRefreshed.v[0];   // a folded doublet: one field for both channels
    rvec_t exc(r.up.size()), vrho(r.up.size());
    for (size_t g=0; g<exc.size(); g++)
    {
        exc[g]=itsXc->GetExcDensity(r.up[g], r.dn[g]);
        vrho[g]=vUp[g]*r.up[g]+vDn[g]*r.dn[g];
    }
    const double dilation=itsSampler->Integrate(exc)-itsSampler->Integrate(vrho);
    // The orbital side, block by block and spin by spin.
    const rmat3d_t orb=D->DM_Virial(this, st);
    for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) V(a,b)+=orb(a,b);
    for (size_t a=1;a<=3;a++) V(a,a)+=dilation;
}

template std::string tVxc_Quadrature<double>::MissingVirial(const tbs_t<double>*, const Structure*) const;
template std::string tVxc_Quadrature<dcmplx>::MissingVirial(const tbs_t<dcmplx>*, const Structure*) const;
template void tVxc_Quadrature<double>::AddVirial(rmat3d_t&, const tDM_CD<double>*, const Structure*) const;
template void tVxc_Quadrature<dcmplx>::AddVirial(rmat3d_t&, const tDM_CD<dcmplx>*, const Structure*) const;
template rmat3d_t tVxc_Quadrature<double>::VirialBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;
template rmat3d_t tVxc_Quadrature<dcmplx>::VirialBlock<dcmplx>(const tobs_t<dcmplx>*, const Spin&, const hmat_t<dcmplx>&, const Structure*) const;
template rmat3d_t tVxc_Quadrature<dcmplx>::VirialBlock<double>(const tobs_t<double>*, const Spin&, const hmat_t<double>&, const Structure*) const;

} // namespace qchem::Hamiltonian
