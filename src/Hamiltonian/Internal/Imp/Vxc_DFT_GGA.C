// File: Hamiltonian/Internal/Imp/Vxc_DFT_GGA.C  The GGA XC term's bodies (GGA stage 2; the declarations and the physics are in the interface).
module;
#include <cassert>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.Vxc_DFT_GGA;
import qchem.Energy;
import qchem.ChargeDensity;

namespace qchem::Hamiltonian
{

template <class T> tVxc_DFT_GGA<T>::tVxc_DFT_GGA(gga_t xc, sampler_t sampler, SpinGroup g) requires std::is_same_v<T,double>
    : itsXc(std::move(xc))
    , itsSampler(std::move(sampler))
    , itsGroup(g)
{
    if (!itsXc)      throw std::invalid_argument("tVxc_DFT_GGA: no GGA functional");
    if (!itsSampler) throw std::invalid_argument("tVxc_DFT_GGA: no density-gradient sampler");
}

// sigma_{st} = grad rho_s . grad rho_t at every point, beside the channel densities: the functional's whole input, one record.
template <class T> GGAInput tVxc_DFT_GGA<T>::InputAt(const tChargeDensity<T>* cd) const
{
    const rvec_t&     up=itsSampler->RhoPol(cd, Spin::Up);
    const rvec_t&     dn=itsSampler->RhoPol(cd, Spin::Down);
    const rvec3vec_t& gu=itsSampler->GradRhoPol(cd, Spin::Up);
    const rvec3vec_t& gd=itsSampler->GradRhoPol(cd, Spin::Down);
    const size_t n=up.size();
    if (dn.size()!=n || gu.size()!=n || gd.size()!=n)
        throw std::logic_error("tVxc_DFT_GGA: the sampler's channels and gradients are not one entry per point");
    rvec_t uu(n), ud(n), dd(n), u(n), d(n);
    for (size_t g=0; g<n; g++)
    {
        // A negative channel (a seed's ringing) is clamped to 0 by Libxc_GGA itself, the rule Libxc_LDA follows; sigma is passed as given.
        u[g]=up[g]; d[g]=dn[g];
        uu[g]=gu[g]*gu[g]; ud[g]=gu[g]*gd[g]; dd[g]=gd[g]*gd[g];
    }
    return GGAInput(std::move(u), std::move(d), std::move(uu), std::move(ud), std::move(dd));
}

template <class T> void tVxc_DFT_GGA<T>::RefreshForDensity(const tChargeDensity<T>* cd) const
{
    assert(cd);
    itsSampler->WarmForDensity(cd, true);   // always spin-resolved: the sampler's collapse makes an unpolarized run its ζ = 0 case
    RequireGradientSampled(InputAt(cd), "the molecular delta route");   // the forgotten-gradient detector (user, 2026-11-20): live, per refresh
}

template <class T> void tVxc_DFT_GGA<T>::BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const
{
    if (cd && bs)
    {
        const GGAInput in=InputAt(cd);
        GGAOutput out(in);
        itsXc->Evaluate(in, out);
        const rvec3vec_t& gu=itsSampler->GradRhoPol(cd, Spin::Up);
        const rvec3vec_t& gd=itsSampler->GradRhoPol(cd, Spin::Down);
        const size_t n=in.size();
        Refreshed& R=itsRefreshed;
        R.v[0]=out.vUp; R.v[1]=out.vDn;
        R.w[0].resize(n); R.w[1].resize(n);
        for (size_t g=0; g<n; g++)   // w_s = 2 v_sigma_ss grad rho_s + v_sigma_ud grad rho_(other)
        {
            R.w[0][g]=(2.0*out.vsUU[g])*gu[g]+out.vsUD[g]*gd[g];
            R.w[1][g]=(2.0*out.vsDD[g])*gd[g]+out.vsUD[g]*gu[g];
        }
        R.source=itsSampler->RhoSource();
        R.densityIsRhoD=cd->HasNoFitComponent();
        R.id=DensityID(cd);
    }
    tDynamic_HT_Imp<T>::BuildBlocks(bs, cd);
}

// A block READS the refreshed field: Up for a folded (unpolarized) block, as the LDA term reads it -- the zeta = 0 collapse gives both
// channels the same field, and H = dE/dD for the total D is either channel's.
template <class T> hmat_t<T> tVxc_DFT_GGA<T>::MakeMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id) const
{
    const Refreshed& R=itsRefreshed;
    if (id.IsNone() || !(id==R.id))
        throw std::logic_error("tVxc_DFT_GGA: block asked for a density the refresh did not see (RefreshForDensity first)");
    itsFockSource=R.source;
    itsFockDensityIsRhoD=R.densityIsRhoD;
    const int c = s==Spin::Down ? 1 : 0;
    return itsSampler->Matrix(bs, R.v[c], R.w[c]);
}

template <class T> void tVxc_DFT_GGA<T>::GetEnergy(EnergyBreakdown& te, const tDM_CD<T>* cd) const
{
    const GGAInput in=InputAt(cd);
    GGAOutput out(in);
    itsXc->Evaluate(in, out);
    te.AddXC(itsSampler->Integrate(out.e), std::nullopt,
             itsFockSource ? ChargeDensity::CertificateLabel(*itsFockSource, itsFockDensityIsRhoD) : "",
             itsFockSource && ChargeDensity::CertifiesRhoD(*itsFockSource, itsFockDensityIsRhoD));
    rvec_t rho(in.size());
    for (size_t g=0; g<rho.size(); g++) rho[g]=in.up[g]+in.dn[g];
    te.charge.lost=itsSampler->Integrate(rho)-cd->GetTotalCharge();   // the mesh-charge leak, as the LDA term reports it
}

template <class T> std::ostream& tVxc_DFT_GGA<T>::Write(std::ostream& os) const
{
    os << "    XC-mesh GGA e(rho_up,rho_dn,sigma) "
       << (itsGroup==SpinGroup::Polarized ? "SPIN-NATIVE" : "(unpolarized: the zeta=0 collapse of the spin-native call)")
       << ", v and w in ONE gather (" << itsSampler->NumPoints() << " points): ";
    itsXc->Write(os);
    return os << std::endl;
}

template class tVxc_DFT_GGA<double>;   // stage 2: the molecular point route

} //namespace
