// File: Hamiltonian/Internal/Imp/PWTerms_XC.C  the XC TERMS -- the thin layer that owns the FUNCTIONAL and nothing else.
//
// One implementation unit of module qchem.Hamiltonian.Internal.PWTerms.  Split 2026-09-08 out of a
// single 1213-line Imp/PWTerms.C (user: "PWTerms.C is huge, again doing too many things") into the
// interface-plus-many-Imp-units shape Internal/Terms.C has always had.  Helpers shared by more than
// one unit (NarrowExact, SampledField) live in the module INTERFACE's non-exported section, which is
// exactly what module-internal linkage is for -- they are visible to every unit of this module and to
// nothing outside it.
module;
#include <mutex>
#include <algorithm>   // std::min (the threaded quadrature's output-column blocking)
#include <cmath>       // std::fabs (the D-MINDRHO residual)
#include <cassert>
#include <complex>
#include <cstdlib>
#include <exception>   // std::exception_ptr (throw containment across the threaded Phi build)
#include <functional>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <memory>
#include <optional>    // the conditionally-charged sub-buckets of the H_xc quadrature
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.Orbital_PP_IBS;             // cast bs ACROSS to the species-field integral service (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.PointSampling;                   // SamplingTarget / SamplingAccuracy (the rho-tilde sampler's tolerance, step 2)


namespace qchem::Hamiltonian
{

// ---- Vxc_Quadrature ------------------------------------------------------------------------------------------

// Built with the SHARED quadrature engine (mesh + Phi tables + per-serial rho) and the imposed spin subgroup.
template <class T> tVxc_Quadrature<T>::tVxc_Quadrature(const xc_t& xc, sampler_t quad, SpinGroup g)
    : itsXc(xc)
    , itsSampler(std::move(quad))
    , itsGroup(g)
{
    assert(itsXc);
    assert(itsSampler);
}

// The periodic term hands the base the same sampler it keeps for the real-TRIM and response faces: ONE engine, two views.
Vxc_Quadrature::Vxc_Quadrature(const xc_t& xc, std::shared_ptr<const ChargeDensity::DensitySampler> quad, SpinGroup g)
    : tVxc_Quadrature<dcmplx>(xc, quad, g)
    , itsPeriodic(std::move(quad))
{}

// THE CHANNEL RASTERS -- the one place the imposed subgroup is read on the Fock/energy path.  Polarized:
// the sampler's {up,dn} pair (a spin-agnostic seed collapses to rho/2 per channel INSIDE RhoPol).  SU(2):
// the folded doublet's one raster, halved into `scratch` and handed over as BOTH channels -- the exact
// zeta=0 collapse, so the functional's spin-native face reproduces the old scalar path bit for bit
// (2*(rho/2)==rho; VWN routes exact zeta=0 through its scalar formula).
template <class T> typename tVxc_Quadrature<T>::Rasters tVxc_Quadrature<T>::ChannelRasters(const tChargeDensity<T>* cd, rvec_t& scratch) const
{
    if (itsGroup==SpinGroup::Polarized)
        return {itsSampler->RhoPol(cd,Spin::Up), itsSampler->RhoPol(cd,Spin::Down)};
    scratch=0.5*itsSampler->Rho(cd);
    return {scratch, scratch};
}

// v^sigma(rho_up, rho_dn) pointwise on the channel rasters, then the engine's Phi-table quadrature (one
// GEMM).  ONE body for both block scalars (Step 3c): the rasters are block-independent; only the final
// quadrature (typed Phi table) differs.  A folded-doublet block (Spin::None) is the same field either
// channel would see, so it asks for Up by convention.
template <class T> rvec_t tVxc_Quadrature<T>::PotentialAt(const Spin& fs, const Rasters& r) const
{
    rvec_t v(r.up.size());
    for (size_t g=0; g<v.size(); g++) v[g]=itsXc->GetVxc(r.up[g], r.dn[g], fs);
    return v;
}
// The block path (V1.41 stage 2, ruling C): it READS what the refresh built for \a cd and computes nothing from rho.  A
// miss is a caller that skipped the refresh, a broken contract no caller can act on, so it throws.
template <class T> const typename tVxc_Quadrature<T>::Refreshed& tVxc_Quadrature<T>::RefreshedFor(DensityID id) const
{
    if (id.IsNone() || !(id==itsRefreshed.id))
        throw std::logic_error("tVxc_Quadrature: block asked for a density the refresh did not see (RefreshForDensity first)");
    // What THIS V_xc saw (pin 28): recorded by the FOCK pass, which is the only pass that asks for blocks, so the energy
    // pass at rho_out reports the source and group residual of the potential the Fock was built from.
    std::lock_guard<std::mutex> lock(itsFockSourceLock.m);
    itsFockSource=itsRefreshed.source;
    itsFockGroupResidual=itsRefreshed.groupResidual;
    return itsRefreshed;
}
template <class T> hmat_t<T> tVxc_Quadrature<T>::MakeMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id) const
{
    const Refreshed& R=RefreshedFor(id);
    const int c = s==Spin::Down ? 1 : 0;
    SizeFourierSampling(s, R, [&](const rvec_t& f){return itsSampler->Matrix(bs, f);});
    return itsSampler->Matrix(bs, R.v[c]);
}
rsmat_t Vxc_Quadrature::MakeMatrix(const robs_t* bs, const Spin& s, DensityID id) const
{
    const Refreshed& R=RefreshedFor(id);
    const int c = s==Spin::Down ? 1 : 0;
    SizeFourierSampling(s, R, [&](const rvec_t& f){return itsPeriodic->Matrix(bs, f);});
    return itsPeriodic->Matrix(bs, R.v[c]);
}

// THE ONE TOLERANCE tau on max|delta rho| of a Fourier series sampled on my points (the rho-tilde sampler, step 2; doc/Pins.md pin 31).
// RIGOROUS, per consumer of the sampled field, each prefactor an integral this quadrature answers:
//   - the XC energy:     |dE_xc| <= tau Sum_sigma Integral |v_xc^sigma|                 => tau_E = target / P_E;
//   - a matrix element:  |dH_ij^sigma| <= tau max_i <i| Sum_t |f_xc^{sigma t}| |i>       => tau_H = target / P_H
//     (Cauchy-Schwarz on a non-negative weight: |<i|w|j>| <= sqrt(<i|w|i><j|w|j>) <= max_i <i|w|i>);
//   - the SCF residual:  the L1 distance of a sampled pair moves by <= 2 tau Integral 1, per electron of the smaller channel, which
//     may use the share kappa_R of the density gate (a CALIBRATION, GridAccuracy::samplerResidualMargin) => tau_R.
// An energy run takes min(tau_E, tau_R), a geometry run min(tau_E, tau_H, tau_R).  ON THE SEED, once (user ruling 2026-11-04: never a
// per-iteration input): every block of the first Fock build lowers it, and the sampler's first SIZED Load freezes it.
template <class T> template <class Mat>
void tVxc_Quadrature<T>::SizeFourierSampling(const Spin& s, const Refreshed& R, Mat matrixOf) const
{
    const rvec_t& v=R.v[s==Spin::Down ? 1 : 0];
    SamplingTarget* tgt=itsSampler->FourierTarget();
    if (!tgt || tgt->loads>0 || !(tgt->accuracy.target>0.0)) return;   // no stated accuracy, or tau already frozen by its first use
    std::lock_guard<std::mutex> lock(itsSamplingLock.m);
    const SamplingAccuracy& A=tgt->accuracy;
    const bool pol=itsGroup==SpinGroup::Polarized;
    rvec_t av(v.size());
    for (size_t g=0; g<v.size(); g++) av[g]=std::fabs(v[g]);
    itsSeedPE[s==Spin::Down ? 1 : 0]=itsSampler->Integrate(av);
    const double PE = pol ? itsSeedPE[0]+itsSeedPE[1] : itsSeedPE[0];
    double tau=A.target/PE;
    std::string from="the XC energy";
    if (A.perElement)
    {
        const auto M=matrixOf(R.seedFxc[s==Spin::Down ? 1 : 0]);
        for (size_t i=0; i<M.rows(); i++) itsSeedPH=std::max(itsSeedPH, std::abs(M(i,i)));
        if (itsSeedPH>0.0 && A.target/itsSeedPH<tau) {tau=A.target/itsSeedPH; from="the XC matrix elements";}
    }
    if (A.residualGate>0.0)
    {
        const double N=R.seedN;
        const double W=itsSampler->Integrate(rvec_t(v.size(), 1.0));
        if (N>0.0 && W>0.0)
            if (const double tR=A.residualMargin*A.residualGate*N/(2.0*W); tR<tau) {tau=tR; from="the SCF residual";}
    }
    // The consumer's target, and the SIZING tolerance: the run's ONE tightening (decided at an SCF boundary by SolidCalculation, never here)
    // shrinks the latter only.
    if (!(tgt->tauTarget>0.0) || tau<tgt->tauTarget) {tgt->tauTarget=tau; tgt->tau=tau*A.tauScale; tgt->tauFrom=from;}
}

template <class T> void tVxc_Quadrature<T>::ReportFourierSampling() const
{
    const SamplingTarget* tgt=itsSampler->FourierTarget();
    if (!tgt || tgt->loads==0 || itsSamplingReported) return;
    std::lock_guard<std::mutex> lock(itsSamplingLock.m);
    if (itsSamplingReported) return;
    // ONCE, when the sampler's first sized use froze tau: the settings the whole SCF runs under.  What it DELIVERS is the samp/τ column.
    std::ostringstream os;   // its own stream: the console's ambient precision would round tau away
    os<<"[XC sampling] rho-tilde NUFFT: tau "<<tgt->tauTarget<<" (set by "<<tgt->tauFrom<<" on the seed: P_E "<<itsSeedPE[0]+itsSeedPE[1];
    if (itsSeedPH>0.0) os<<", P_H "<<itsSeedPH;
    os<<"), eta "<<tgt->accuracy.eta;
    if (tgt->accuracy.tauScale<1.0) os<<", sized at the TIGHTENED "<<tgt->tau<<" (x"<<tgt->accuracy.tauScale<<")";
    std::cout<<os.str()<<std::endl;
    qchem::report::Scope here("xcSampling");
    qchem::report::Set("tau", tgt->tauTarget);
    qchem::report::Set("tauSizing", tgt->tau);
    qchem::report::Set("tauFrom", tgt->tauFrom);
    qchem::report::Set("seedPE", itsSeedPE[0]+itsSeedPE[1]);
    if (itsSeedPH>0.0) qchem::report::Set("seedPH", itsSeedPH);
    qchem::report::Set("eta", tgt->accuracy.eta);
    qchem::report::Set("tauScale", tgt->accuracy.tauScale);
    itsSamplingReported=true;
}

template <class T> void tVxc_Quadrature<T>::GetEnergy(EnergyBreakdown& te, const tDM_CD<T>* cd) const
{
    rvec_t half;
    const Rasters r=ChannelRasters(cd, half);   // reuses the iteration's rasters (same density serial)
    rvec_t exc(r.up.size()), rho(r.up.size());
    // E_xc = Integral e(rho_up,rho_dn), the PER-VOLUME energy density, because that is the form that
    // COMPOSES: exchange contributes Sum_sigma eps_x(rho_sigma) rho_sigma and correlation eps_c rho_tot, and
    // those share no denominator.
    for (size_t g=0; g<exc.size(); g++) { exc[g]=itsXc->GetExcDensity(r.up[g], r.dn[g]); rho[g]=r.up[g]+r.dn[g]; }
    te.Add("Exc", itsSampler->Integrate(exc), EnergyRole::Potential);
    // The charge accounting of the same density: the mesh-charge leak (the quadrature's health metric --
    // CP2K's grid-charge-lost readout) is the quadrature integral of rho vs the analytic Tr(DS); the
    // integrated site moments are the observable, where both rasters are in hand (R1.0h).
    te.charge.lost        = itsSampler->Integrate(rho) - cd->GetTotalCharge();
    te.charge.siteMoments = SiteMoments(cd);
    if (itsFockSource)
    {
        te.charge.xcRho       = ChargeDensity::ShortName(*itsFockSource);
        te.charge.xcRhoFitted = *itsFockSource==ChargeDensity::XCRhoSource::Fitted;
    }
    te.charge.xcGroupResidual = itsFockGroupResidual;
    ReportFourierSampling();
    // THE samp/τ COLUMN: the worst audited delivered/tau of the sized samplings since the last energy that saw new ones (one iteration's:
    // the Fock build's correction and the residual pair).  An energy pass with no new samplings (GDM's line search) repeats the value.
    if (SamplingTarget* tgt=itsSampler->FourierTarget(); tgt && tgt->loads>0)
    {
        std::lock_guard<std::mutex> lock(itsSamplingLock.m);
        if (tgt->loads!=itsSamplingLastLoads) {itsSamplingShown=tgt->TakeWindow(); itsSamplingLastLoads=tgt->loads;}
        te.charge.samplingOverTau=itsSamplingShown;
    }
}

// THE OBSERVABLE (doc/OpenWork.md N1/T2; R1.0h): mu_A = Integral w_A (rho_up - rho_dn), in electrons, over the
// sampler's atom-centred partition.  FREE -- both channel rasters are cached for this density serial -- and
// the term's, not the sampler's, because the term is what knows the field is a spin difference; the
// sampler only knows how to partition-integrate a field at its points.  Empty when the run resolves no
// spin (m == 0 identically) or the quadrature has no site blocks -- the caller treats empty as "this run
// cannot answer", never as "zero".
template <class T> rvec_t tVxc_Quadrature<T>::SiteMoments(const tChargeDensity<T>* cd) const
{
    if (itsGroup!=SpinGroup::Polarized) return rvec_t();
    return itsSampler->SiteIntegrals(rvec_t(itsSampler->RhoPol(cd,Spin::Up)-itsSampler->RhoPol(cd,Spin::Down)));
}

// D-MINDRHO: the sampler caches ONE density serial (both channels of it), so BOTH of a's channels are copied before b's
// replace them -- two rasters per call, not four (the a-up, b-up, a-dn, b-dn order re-rastered each density).
template <class T> double tVxc_Quadrature<T>::L1Distance(const tChargeDensity<T>* a, const tChargeDensity<T>* b) const
{
    auto l1=[&](const rvec_t& fa, const rvec_t& fb)
    {
        rvec_t d(fa.size());
        for (size_t g=0; g<d.size(); g++) d[g]=std::fabs(fa[g]-fb[g]);
        return itsSampler->Integrate(d);
    };
    if (itsGroup!=SpinGroup::Polarized)
    {
        const rvec_t fa=itsSampler->Rho(a);
        return l1(fa, itsSampler->Rho(b));
    }
    const rvec_t au=itsSampler->RhoPol(a,Spin::Up), ad=itsSampler->RhoPol(a,Spin::Down);
    const rvec_t bu=itsSampler->RhoPol(b,Spin::Up), bd=itsSampler->RhoPol(b,Spin::Down);
    return l1(au,bu)+l1(ad,bd);
}

// THE RESPONSE (R2): the ALDA kernel at the linearisation density, times δρ_σ, gathered by the same adjoint.
void Vxc_Quadrature::RefreshForDensity(const cbs_t*, const cChargeDensity* D0, const TransitionDensity<dcmplx>& delta) const
{
    if (!D0) throw std::logic_error("Vxc_Quadrature: the XC response needs the density it linearises about (D0)");
    if (D0->Version()!=itsKernelVersion)
    {
        qchem::report::Timed timed("response: f_xc at rho0");
        rvec_t half;
        const Rasters r=ChannelRasters(D0, half);
        const size_t n=r.up.size();
        itsFuu.resize(n); itsFud.resize(n); itsFdu.resize(n); itsFdd.resize(n);
        for (size_t g=0; g<n; g++)
        {
            itsFuu[g]=itsXc->GetFxc(r.up[g], r.dn[g], Spin::Up,   Spin::Up  );
            itsFud[g]=itsXc->GetFxc(r.up[g], r.dn[g], Spin::Up,   Spin::Down);
            itsFdu[g]=itsXc->GetFxc(r.up[g], r.dn[g], Spin::Down, Spin::Up  );
            itsFdd[g]=itsXc->GetFxc(r.up[g], r.dn[g], Spin::Down, Spin::Down);
        }
        itsKernelVersion=D0->Version();
        itsDeltaVersion=size_t(-1);   // the kernel moved: every δv is stale
    }
    if (delta.Version()==itsDeltaVersion) return;
    cvec_t du, dd;
    if (itsGroup==SpinGroup::Polarized) {du=itsPeriodic->Sample(delta, Spin::Up); dd=itsPeriodic->Sample(delta, Spin::Down);}
    else                                {du=0.5*itsPeriodic->Sample(delta, Spin::None); dd=du;}
    itsDeltaQ=qchem::ChargeDensity::WaveVectorOf(*delta.Coupling());
    const size_t n=du.size();
    if (n!=itsFuu.size()) throw std::logic_error("Vxc_Quadrature: δρ and the kernel live on different quadratures");
    itsDvUp.resize(n); itsDvDn.resize(n);
    for (size_t g=0; g<n; g++)
    {
        itsDvUp[g]=itsFuu[g]*du[g]+itsFud[g]*dd[g];
        itsDvDn[g]=itsFdu[g]*du[g]+itsFdd[g]*dd[g];
    }
    itsDeltaVersion=delta.Version();
}

mat_t<dcmplx> Vxc_Quadrature::GetMatrix(const cobs_t* bra, const cobs_t* ket, const Spin& s, const TransitionDensity<dcmplx>& delta) const
{
    if (delta.Version()!=itsDeltaVersion)
        throw std::logic_error("Vxc_Quadrature: GetMatrix for a transition density the response phase did not refresh");
    // A folded-doublet block (Spin::None) sees the same δv either channel would: Up by convention, as MakeMatrixT.
    return itsPeriodic->Matrix(bra, ket, s==Spin::Down ? itsDvDn : itsDvUp, itsDeltaQ);
}

template <class T> std::ostream& tVxc_Quadrature<T>::Write(std::ostream& os) const
{
    return os << "    XC-mesh " << (itsGroup==SpinGroup::Polarized ? "SPIN-NATIVE v_xc^sigma(rho_up,rho_down)"
                                                                   : "v_xc(rho(r))")
              << ", exchange+correlation in ONE gather (" << itsSampler->NumPoints() << " points)." << std::endl;
}

// ---- THE EAGER REFRESH PHASE (doc/OpenWork.md item KP) -------------------------------------------------
// rho on the quadrature's points is k-INDEPENDENT -- one array (or one pair), correct for every Bloch block
// -- but it was filled by whichever block asked first.  Warming it here makes the block loop read-only in
// the ordinary path, which is the precondition for running the blocks concurrently.
// ★ V1.41 STAGE 2: the refresh also builds v^sigma, once per spin channel, which every block then only reads (it used to be
// recomputed per block from the rasters).  And on the seed, while tau is unfrozen, the |f_xc| fields and the charge the sizing needs.
// Every refresh, the energy pass's included: warm the channel raster(s) GetEnergy integrates.  The potential is not built here:
// only the blocks read it, so BuildPotential builds it in the Fock refresh's BuildBlocks (the energy pass, bs = null, never pays
// for v^sigma on the whole mesh: +2.4 s of a 247 s MnO run when it did, lane S 2026-11-07).
template <class T> void tVxc_Quadrature<T>::RefreshForDensity(const tChargeDensity<T>* cd) const
{
    assert(cd);
    itsSampler->WarmForDensity(cd, itsGroup==SpinGroup::Polarized);
}

template <class T> void tVxc_Quadrature<T>::BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const
{
    BuildPotential(bs, cd);
    tDynamic_HT_Imp<T>::BuildBlocks(bs, cd);
}

// The block-only state for cd, built at the Fock refresh (RefreshForDensity(cd) has warmed the rasters just before, so
// ChannelRasters reads the sampler's cache).  bs null: an energy evaluation, no blocks, nothing to build.
template <class T> void tVxc_Quadrature<T>::BuildPotential(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const
{
    if (!cd || !bs) return;
    const bool pol=itsGroup==SpinGroup::Polarized;
    rvec_t half;
    const Rasters r=ChannelRasters(cd, half);
    Refreshed& R=itsRefreshed;
    R.source=itsSampler->RhoSource();          // read while the sampler's cache is cd's
    R.groupResidual=itsSampler->GroupResidual();
    const int nc = pol ? 2 : 1;               // a folded doublet (Spin::None) is the same field either channel would see: Up
    for (int c=0; c<nc; c++) R.v[c]=PotentialAt(c==1 ? Spin::Down : Spin::Up, r);
    for (int c=nc; c<2; c++) R.v[c]=rvec_t();
    const SamplingTarget* tgt=itsSampler->FourierTarget();
    const bool sizing = tgt && tgt->loads==0 && tgt->accuracy.target>0.0;
    for (int c=0; c<2; c++) R.seedFxc[c]=rvec_t();
    R.seedN=0.0;
    if (sizing && tgt->accuracy.perElement)
        for (int c=0; c<nc; c++)
        {
            const Spin fs = c==1 ? Spin::Down : Spin::Up;
            rvec_t& f=R.seedFxc[c];
            f.resize(r.up.size());
            for (size_t g=0; g<f.size(); g++)
                f[g] = pol ? std::fabs(itsXc->GetFxc(r.up[g], r.dn[g], fs, Spin::Up))+std::fabs(itsXc->GetFxc(r.up[g], r.dn[g], fs, Spin::Down))
                           : 0.5*std::fabs(itsXc->GetFxc(r.up[g], r.dn[g], Spin::Up, Spin::Up)+itsXc->GetFxc(r.up[g], r.dn[g], Spin::Up, Spin::Down));
        }
    if (sizing && tgt->accuracy.residualGate>0.0)
    {
        double N=cd->GetTotalCharge();
        if (pol)
            if (const auto* up=ChargeDensity::ChannelOf(cd,Spin::Up); up)
                if (const auto* dn=ChargeDensity::ChannelOf(cd,Spin::Down); dn) N=std::min(up->GetTotalCharge(), dn->GetTotalCharge());
        R.seedN=N;
    }
    R.id=DensityID(cd);
}

template class tVxc_Quadrature<double>;   // the molecular / atom path (V-XC-ONE-TERM, 2026-10-31)
template class tVxc_Quadrature<dcmplx>;   // the periodic path (Vxc_Quadrature)

} //namespace
