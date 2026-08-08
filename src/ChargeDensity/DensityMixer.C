// File: DensityMixer.C  The DENSITY-FACE of SCF convergence (doc/SCFStrategyPlan.md).
//
// A density mixer folds a freshly diagonalised density (rho_out) into the running one (rho_in) and returns
// the convergence gate ‖Δρ‖.  It is one of the four SCF role-seams; the SCFIterator owns the density
// LIFECYCLE (SetWorkingCD / lineage / GetTotalEnergy) and the mixer owns the POLICY + arithmetic + state.
//
// Concretes:
//   * LinearMixer(α)  -- ρ_next = (1−α)ρ_in + α ρ_out via IrrepCD::MixIn, plus the historical [F,D]-keyed
//                        adaptive-α re-damp/grow.  α=1 is passthrough, so there is NO NullMixer (§3 of the
//                        plan): "no mixing" == the molecular default (StartingRelaxRo defaults to 1.0).
//   * KerkerMixer(α,G0) -- Kerker-preconditioned ρ̃(G)-mixing on FourierMixCD (periodic / dcmplx only).
//
// This is Increment 1 of the plan: a BEHAVIOUR-PRESERVING extraction of the mixing that was inlined in
// tSCFIterator::Iterate.  The face is deliberately shaped to reproduce the legacy control flow bit-for-bit
// (the re-damp interleaves with the energy evaluation, which the iterator still drives); it will slim down
// when the shared extrapolator + Pulay land and the adaptive-α wart is reconsidered.
module;
#include <memory>
#include <iostream>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <deque>
#include <type_traits>
#include <optional>
#include <cstdlib>   // getenv (the QCHEM_SPINBLIND_KERKER A/B valve)
export module qchem.ChargeDensity.DensityMixer;
export import qchem.ChargeDensity;                 // tChargeDensity<T>, tDM_CD<T>
import qchem.ChargeDensity.FourierDensity;         // FourierDensity, ΔG_Map
import qchem.ChargeDensity.FourierMixCD;           // FourierMixCD / KerkerMix
import qchem.Math.DIIS;                             // the shared Pulay/DIIS bordered-solve engine
import qchem.Blaze;                                 // rsmat_t/rvec_t/ivec3_t + blazem::zero (the Pulay B-solve)
import qchem.BasisSet;                             // tBasisSet<T>, operator[]
import qchem.BasisSet.Band_FT_IBS;                 // Band_FT_IBS::CreateVxcFitBasisSet
import qchem.BasisSet.Fit_IBS;                     // cFIT_SF_ABS (the FourierDensity face arg)
import qchem.BasisSet.G_FieldEvaluator;            // ApplySpectralFilter (the raster Kerker step, 0.5(f2))
import qchem.ReciprocalLattice;                    // ReciprocalLattice
import qchem.UnitCell;                             // UnitCell (+ MakeReciprocalCell)
import qchem.Structure;                            // Structure
import qchem.Mesh;                                 // qcMesh::MeshParams
import qchem.Types;                                // dcmplx

export namespace qchem::ChargeDensity
{

//! Signals available to an adaptive density-mixing policy once the energy + orbital gradient are known.
struct MixSignals { double E=0.0, FD=0.0, FDold=0.0; };

//! The density-face of SCF convergence.  The SCFIterator calls, per fixed-point iteration:
//!   working = fresh diagonalised density (already made the lineage head by the iterator)
//!   dρ = Mix(working, old);            // fold rho_out into the running density
//!   ... iterator computes E, [F,D] ...
//!   if (WantsReDamp(sig)) { iterator reseats working->fresh; dρ = ReDampMix(working, old); recompute E; }
//!   UpdateRelax(sig);
//! and drives the next Fock from FockDensity(working).  Non-adaptive mixers (Kerker, future Pulay) take the
//! no-op defaults for the three adaptive hooks.
template <class T> class tDensityMixer
{
public:
    typedef std::shared_ptr<tDM_CD<T>> cd_t;
    virtual ~tDensityMixer() {}
    //! Fold the fresh \a working density into the running one; returns ‖Δρ‖ (the convergence gate).
    virtual double Mix(cd_t& working, const cd_t& old) = 0;
    //! The density that drives the NEXT Fock (default: the working density; Kerker → the mixed ρ̃).
    virtual const tChargeDensity<T>* FockDensity(const cd_t& working) const { return working.get(); }
    //! The current step size α (for the SCF trace only).
    virtual double GetRelax() const = 0;
    //! A 3-char self-identifier for the per-iteration ρ_mix column (doc/GPWPlan1.md item 2): "Lin"
    //! (linear D-mixing), "Ker" (Kerker ρ̃), "Pul" (density-DIIS/Pulay).  Default = the linear mixer.
    virtual const char* Tag() const { return "Lin"; }
    //! Adaptive [F,D]-keyed policy (LinearMixer only; no-op elsewhere).  Post-energy re-damp on divergence.
    virtual bool   WantsReDamp(const MixSignals&) const { return false; }
    //! Re-mix \a working (already reseated to the fresh density by the iterator) more aggressively; ‖Δρ‖.
    virtual double ReDampMix(cd_t& /*working*/, const cd_t& /*old*/) { return 0.0; }
    //! Grow/clamp the step for the next iteration.
    virtual void   UpdateRelax(const MixSignals&) {}
};

//! Linear density-matrix mixing with the legacy adaptive-α heuristics.  α=1 (the default) = passthrough.
template <class T> class LinearMixer : public tDensityMixer<T>
{
public:
    typedef typename tDensityMixer<T>::cd_t cd_t;
    explicit LinearMixer(double relax0) : itsRelax(relax0) {}

    double Mix(cd_t& working, const cd_t& old) override
    {
        double dcd = working->GetChangeFrom(*old)/working->GetTotalCharge();   // relative MaxAbs change
        if (dcd<1e-5) itsRelMax=0.5;
        working->MixIn(*old, 1.0-itsRelax);                                    // (1−relax)ρ_in + relax ρ_out
        return dcd;
    }
    double GetRelax() const override { return itsRelax; }

    bool WantsReDamp(const MixSignals& s) const override
    {
        double dFD = s.FD - s.FDold;
        return s.FD>s.FDold && std::fabs(dFD)>1e-9;
    }
    double ReDampMix(cd_t& working, const cd_t& old) override
    {
        double dcd = working->GetChangeFrom(*old);      // NOTE: un-normalised -- matches the legacy re-damp
        working->MixIn(*old, 1.0-itsRelax/4.0);
        itsRelax*=0.8;
        return dcd;
    }
    void UpdateRelax(const MixSignals& s) override
    {
        if (s.FD<s.FDold) itsRelax*=1.5;
        if (itsRelax>itsRelMax) itsRelax=itsRelMax;
    }
private:
    double itsRelax;
    double itsRelMax=1.0;
};

// --- The RASTER-SPACE Kerker step (the raw-XC shadow of FourierMixCD::KerkerMix; doc/GPWPlan 0.5(f2)) ---
// rho_mix(r) = rho_in + alpha * F^-1[ G^2/(G^2+G0^2) F(rho_out - rho_in) ] over the FULL box.  The filter is
// a SMOOTH multiplier (no truncation -> no Gibbs), k(0)=0 conserves charge -- identical algebra to the ball
// mix but over every mode the raster represents, so the raw feed's out-of-ball content keeps mixing at k~1.
inline rvec_t RasterKerker(const BasisSet::G_FieldEvaluator& ge, const rvec_t& in, const rvec_t& out,
                           double alpha, double G0)
{
    const double G0sq=G0*G0;
    rvec_t delta=out; delta-=in;
    rvec_t mix=in;
    mix+=alpha*ge.ApplySpectralFilter(delta, [G0sq](double g2){return g2/(g2+G0sq);});
    return mix;
}

//! THE MIXABLE ENTITY: a G-space field \f$\tilde f(G)\f$ plus its optional raw-raster shadow.
//!
//! This type exists to make \f$m=\rho_\uparrow-\rho_\downarrow\f$ a FIRST-CLASS SUBJECT OF MIXING (user,
//! 2026-08-08).  A ρ̃ mixer's real precondition was never "a density" -- Kerker and Pulay only ever touched
//! the \c FourierDensity face -- yet they DECLARED \c Mix(cd_t&,...), i.e. a density-matrix-backed density.
//! That wider-than-necessary precondition is precisely what made m unsubstitutable: m has no density matrix,
//! no charge and no positivity, so it could never satisfy a requirement none of those mixers actually had.
//! Naming the true subject fixes it -- ρ and m are now THE SAME TYPE, and any field mixer takes either.
//! (Faking a density instead -- a \c tDM_CD with asserting stubs -- is the anti-pattern this codebase
//! avoids: capabilities live only on the types that have them.)
//!
//! NB the asymmetry is real and deliberate: the mix INPUT is a field (ρ or m equally), while the Fock
//! OUTPUT is genuinely a density -- the Hamiltonian consumes ρ, never m.  So substitutability belongs on
//! the input side, which is where it now lives.
struct GField
{
    ΔG_Map tilde;    //!< \f$\tilde f(G)\f$ on the fit basis's difference set
    rvec_t raster;   //!< the raw-raster shadow (doc/GPWPlan 0.5(f2)); empty = the raw pipeline is off
};

//! A mixer whose subject is a \c GField -- i.e. one that works purely in G space (Kerker, Pulay; NOT the
//! DM-based linear mixer, which genuinely needs the density matrix).  \c tDensityMixer::Mix is then just the
//! adapter that extracts the field from a density and calls this.
class tFieldMixer
{
public:
    virtual ~tFieldMixer() {}
    //! Fold the freshly collocated \a out into the running mixed field; returns the residual ‖out−in‖_∞.
    virtual double MixField(const GField& out) = 0;
    //! The running mixed density -- the PRESENTATION of the mixed field (what \c FockDensity returns, and
    //! what a caller recombining spin channels reads ρ̃ back from).
    virtual const FourierMixCD& Mixed() const = 0;
};

//! Kerker-preconditioned ρ̃(G)-mixing (periodic / dcmplx).  ρ_mix = ρ_in + α·G²/(G²+G0²)·(ρ_out−ρ_in).
//! NB \b G0=0 makes this the PLAIN LINEAR G-space mixer (the filter is identically 1) -- which is exactly
//! the "linear on m" leaf of the (ρ,m) channel basis, so that construction needs no new mixer type.
//! Holds the running mixed ρ̃ as a FourierMixCD; the next Fock is driven from it.  Built by
//! MakeGSpaceMixer -- on a polarized run, one of these PER SPIN CHANNEL (see PolarizedDensityMixer).
//! NB G=0 IS mixed, at full α: our ρ̃ is a fit-basis PROJECTION whose (0,0,0) coefficient is shape-dependent
//! rather than the fixed N/Ω, so freezing it would strand the XC's mean density at the seed (the reason is
//! in FourierMixCD::KerkerMix, which owns the filter).  CP2K does the same -- its kerker_factor array is
//! left at 1.0 for G=0 (qs_mixing_utils.F: `ig1 = 2` when the grid has G=0).
class KerkerMixer : public tDensityMixer<dcmplx>, public virtual tFieldMixer
{
public:
    KerkerMixer(double relax, double G0, std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
                std::shared_ptr<FourierMixCD> rho0, rvec_t raw0)
        : itsRelax(relax), itsKerkerG0(G0), itsKerkerFit(std::move(fit)), itsMixedRho(std::move(rho0))
        , itsRawIn(std::move(raw0))
    { if (itsRawIn.size()) itsMixedRho->SetRawRho(itsRawIn); }

    //! The density-face entry: extract ρ̃_out + the raster shadow, then do the G-space arithmetic below.
    double Mix(cd_t& working, const cd_t& /*old*/) override
    {
        auto* fd = dynamic_cast<const FourierDensity*>(working.get());
        assert(fd && itsMixedRho);
        return MixField({fd->GetFourierDensity(*itsKerkerFit), fd->GetRhoOnGrid(*itsKerkerFit)});
    }
    const FourierMixCD& Mixed() const override { assert(itsMixedRho); return *itsMixedRho; }

    double MixField(const GField& out) override
    {
        assert(itsMixedRho);
        const ΔG_Map& rho_out = out.tilde;
        const rvec_t& rawOut  = out.raster;
        const ΔG_Map& rho_in  = itsMixedRho->RhoTilde();
        // SCF residual ‖ρ̃_out − ρ̃_in‖_∞ -- the RIGHT ρ-mixing gate (0 at the fixed point).
        double resid = 0.0;
        for (const auto& [dm, ro] : rho_out)
        {
            auto it = rho_in.find(dm);
            resid = std::max(resid, std::abs(dcmplx(ro) - (it!=rho_in.end() ? dcmplx(it->second) : dcmplx(0.0))));
        }
        for (const auto& [dm, ri] : rho_in)
            if (rho_out.find(dm)==rho_out.end()) resid = std::max(resid, std::abs(dcmplx(ri)));
        itsMixedRho.reset(FourierMixCD::KerkerMix(*itsMixedRho, rho_out, itsRelax, itsKerkerG0));
        // RAW-raster shadow (0.5(f2)): the same Kerker step on rho_raw(r), deposited so the XC feed stays raw
        // through the DYNAMICS.  Late-activates the first time the working density answers raw (a SAD-seeded
        // run's iteration 1); deactivates for the run if a raw answer stops coming or changes raster.
        auto*  ge     = dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsKerkerFit.get());
        if (rawOut.size() && ge)
        {
            if (itsRawIn.size()!=rawOut.size()) itsRawIn=rawOut;                       // bootstrap/late-activate
            itsRawIn=RasterKerker(*ge, itsRawIn, rawOut, itsRelax, itsKerkerG0);
            itsMixedRho->SetRawRho(itsRawIn);
        }
        else itsRawIn=rvec_t{};
        return resid;
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t&) const override { return itsMixedRho.get(); }
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "Ker"; }
private:
    double itsRelax, itsKerkerG0;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsKerkerFit;
    std::shared_ptr<FourierMixCD>                itsMixedRho;
    rvec_t itsRawIn;   //!< the raster shadow of itsMixedRho (empty = raw pipeline off)
};

// --- ΔG_Map arithmetic for Pulay (keys = integer G-index, consistent across iterations) ---
inline ΔG_Map MapSub(const ΔG_Map& a, const ΔG_Map& b)   // a - b
{
    ΔG_Map r=a;
    for (const auto& [k,v]:b) r[k]-=v;
    return r;
}
inline ΔG_Map MapCombine(const std::deque<ΔG_Map>& maps, const rvec_t& c)  // Σ cᵢ mapᵢ
{
    ΔG_Map r;
    for (size_t i=0;i<maps.size();++i)
        for (const auto& [k,v]:maps[i]) r[k]+=c[i]*v;
    return r;
}
inline ΔG_Map MapAdd(ΔG_Map a, const ΔG_Map& b)          // a + b (the ↑+↓ channel sum)
{
    for (const auto& [k,v]:b) a[k]+=v;
    return a;
}
inline ΔG_Map MapScale(ΔG_Map a, double s)               // s*a (the ½ of the (ρ±m)/2 channel rebuild)
{
    for (auto& [k,v]:a) v=s*v;
    return a;
}
//! s*(a + sign*b) for the raster shadows, empty unless BOTH arms answer on the same raster.
inline rvec_t RawCombine(const rvec_t& a, const rvec_t& b, double sign, double s)
{
    if (a.size()==0 || a.size()!=b.size()) return rvec_t{};
    rvec_t r=a; r+=sign*b; r*=s;
    return r;
}
inline double MapMaxAbs(const ΔG_Map& m)
{
    double x=0.0; for (const auto& [k,v]:m) x=std::max(x,std::abs(v)); return x;
}
inline double MapInnerRe(const ΔG_Map& a, const ΔG_Map& b)  // Re Σ_{G≠0} conj(aᵢ)·bⱼ (G=0 is never mixed)
{
    double s=0.0;
    for (const auto& [k,v]:a)
    {
        if (k==ivec3_t(0,0,0)) continue;
        auto it=b.find(k);
        if (it!=b.end()) s+=std::real(std::conj(v)*it->second);
    }
    return s;
}

//! PULAY (density-DIIS) ρ̃-mixing, Kerker-preconditioned (periodic / dcmplx).  Keeps a history of the fed
//! densities ρ̃_in and the freshly collocated ρ̃_out; each step solves the DIIS bordered system (the shared
//! qchem.Math.DIIS engine) over the residuals ρ̃_out−ρ̃_in for the optimal coefficients c (Σc=1), forms the
//! extrapolated ρ̃_in*=Σcᵢρ̃_inᵢ and ρ̃_out*=Σcᵢρ̃_outᵢ, and applies the Kerker step to THOSE via FourierMixCD::
//! KerkerMix (= ρ̃_in* + α·G²/(G²+G0²)·(ρ̃_out*−ρ̃_in*)).  First iteration (history<2) falls back to plain
//! Kerker.  doc/SCFStrategyPlan.md §4 (the density-face use of the shared extrapolator).
class PulayMixer : public tDensityMixer<dcmplx>, public virtual tFieldMixer
{
public:
    PulayMixer(double relax, double G0, int depth, int start, std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
               ReciprocalLattice recip, ΔG_Map rho0, double charge, rvec_t raw0)
        : itsRelax(relax), itsKerkerG0(G0), itsDepth(depth), itsStart(start), itsKerkerFit(std::move(fit))
        , itsRecip(recip), itsCharge(charge)
        , itsMixedRho(std::make_shared<FourierMixCD>(std::move(rho0), recip, charge))
        , itsRawIn(std::move(raw0))
    { if (itsRawIn.size()) itsMixedRho->SetRawRho(itsRawIn); }

    //! The density-face entry: extract ρ̃_out + the raster shadow, then do the G-space arithmetic below.
    double Mix(cd_t& working, const cd_t&) override
    {
        auto* fd = dynamic_cast<const FourierDensity*>(working.get());
        assert(fd && itsMixedRho);
        return MixField({fd->GetFourierDensity(*itsKerkerFit), fd->GetRhoOnGrid(*itsKerkerFit)});
    }
    const FourierMixCD& Mixed() const override { assert(itsMixedRho); return *itsMixedRho; }

    double MixField(const GField& field) override
    {
        assert(itsMixedRho);
        const ΔG_Map& out    = field.tilde;
        const rvec_t& rawOut = field.raster;
        ΔG_Map in  = itsMixedRho->RhoTilde();               // ρ̃_in : the density fed to this iteration's Fock
        ΔG_Map res = MapSub(out,in);                        // residual = ρ̃_out − ρ̃_in
        double resid = MapMaxAbs(res);
        // RAW-raster shadow inputs (0.5(f2)); late-activates like KerkerMixer, drops out if answers stop.
        auto*  ge     = dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsKerkerFit.get());
        const bool raw = rawOut.size() && ge;
        if (raw && itsRawIn.size()!=rawOut.size()) itsRawIn=rawOut;                    // bootstrap/late-activate
        if (!raw) { itsRawIn=rvec_t{}; itsRawIns.clear(); itsRawOuts.clear(); }
        rvec_t rawIn = itsRawIn;                            // the shadow of `in` (captured before the update)

        // PRIME with plain Kerker until we are near the fixed point (history-based mixing is unstable far
        // out).  No history is accumulated during priming, so Pulay starts with clean, linear-regime residuals.
        if (++itsCount<=itsStart)
        {
            itsMixedRho.reset(FourierMixCD::KerkerMix(*itsMixedRho,out,itsRelax,itsKerkerG0));
            if (raw)
            {
                itsRawIn=RasterKerker(*ge, rawIn, rawOut, itsRelax, itsKerkerG0);
                itsMixedRho->SetRawRho(itsRawIn);
            }
            return resid;
        }

        itsIns.push_back(in); itsOuts.push_back(out); itsResiduals.push_back(res);   // grow the history
        if (raw) { itsRawIns.push_back(rawIn); itsRawOuts.push_back(rawOut); }
        while ((int)itsResiduals.size()>itsDepth)                                    // prune to `depth`
            { itsIns.pop_front(); itsOuts.pop_front(); itsResiduals.pop_front(); }
        while (itsRawIns.size()>itsResiduals.size()) { itsRawIns.pop_front(); itsRawOuts.pop_front(); }

        const size_t n=itsResiduals.size();
        ΔG_Map inStar, outStar;
        rvec_t rawInStar=rawIn, rawOutStar=rawOut;         // raw stars default to the plain pair
        if (n<2) { inStar=in; outStar=out; }               // not enough history yet → plain Kerker
        else
        {
            rsmat_t B=blazem::zero<double>(n);              // Bᵢⱼ = ⟨resᵢ,resⱼ⟩ (symmetric)
            for (size_t i=0;i<n;++i)
                for (size_t j=i;j<n;++j) B(i,j)=MapInnerRe(itsResiduals[i],itsResiduals[j]);
            rvec_t c=qchem::Math::DIIS::Coefficients(qchem::Math::DIIS::Bordered(B));
            inStar =MapCombine(itsIns ,c);                 // ρ̃_in*  = Σ cᵢ ρ̃_inᵢ
            outStar=MapCombine(itsOuts,c);                 // ρ̃_out* = Σ cᵢ ρ̃_outᵢ
            // The raw shadow takes the SAME extrapolation coefficients -- but only once its history spans the
            // whole window (a late-activated shadow falls back to the plain pair until the deques align).
            if (raw && itsRawIns.size()==n)
            {
                rawInStar =c[0]*itsRawIns[0];  for (size_t i=1;i<n;++i) rawInStar +=c[i]*itsRawIns[i];
                rawOutStar=c[0]*itsRawOuts[0]; for (size_t i=1;i<n;++i) rawOutStar+=c[i]*itsRawOuts[i];
            }
        }
        FourierMixCD inStarCD(inStar,itsRecip,itsCharge);  // Kerker step on the DIIS-extrapolated pair
        itsMixedRho.reset(FourierMixCD::KerkerMix(inStarCD,outStar,itsRelax,itsKerkerG0));
        if (raw)
        {
            itsRawIn=RasterKerker(*ge, rawInStar, rawOutStar, itsRelax, itsKerkerG0);
            itsMixedRho->SetRawRho(itsRawIn);
        }
        return resid;
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t&) const override { return itsMixedRho.get(); }
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "Pul"; }
private:
    double itsRelax, itsKerkerG0; int itsDepth, itsStart, itsCount=0;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsKerkerFit;
    ReciprocalLattice itsRecip;
    double itsCharge;
    std::shared_ptr<FourierMixCD> itsMixedRho;
    std::deque<ΔG_Map> itsIns, itsOuts, itsResiduals;      // the Pulay history (aligned index-wise)
    rvec_t itsRawIn;                                       //!< the raster shadow of itsMixedRho (empty = off)
    std::deque<rvec_t> itsRawIns, itsRawOuts;              //!< its history (aligned with itsIns while active)
};

//---------------------------------------------------------------------------------------------------------
//  THE POLARIZED ρ̃ MIXER (2026-08-07; doc/SymmetryUpgradePlan.md §7 step 7)
//
//  A ρ̃-space mixer carries ONE FourierMixCD -- the ↑+↓ total, no spin channels -- and drives every Fock
//  from it.  XC_GridEngine::RhoPol then finds no spin face and takes its ρ↑=ρ↓=ρ/2 branch, so v_xc^↑≡v_xc^↓
//  and a POLARIZED run is silently unpolarized from iteration 1 (measured on MnO: a seed staggered at
//  m_stag=0.366 read EXACTLY 0 at iteration 1).  The cure is a COMPOSITION, not a second implementation:
//  one ordinary mixer per spin channel, plus the view below that gives their pair both faces the framework
//  consumes.  Every mixer that works through the FourierDensity face -- Kerker today, Pulay today, whatever
//  lands next -- works here unchanged, because the composite never learns which one it holds.
//
//  BASIS NOTE (checked against CP2K's qs_gspace_mixing.F, 2026-08-07).  CP2K transforms to (ρ_total, m)
//  in its mixing DRIVER and then runs the same per-channel loop with the SAME Kerker factor and α on both.
//  For a LINEAR operator that is algebraically identical to mixing (ρ↑,ρ↓) per channel -- Kerker is linear
//  in the residual, so m_mix = m_in + αK(m_out−m_in) either way -- i.e. THIS composite reproduces CP2K's
//  Kerker exactly, with no proof burden.  The basis only becomes a real choice for the NONLINEAR history
//  mixers (Pulay/Broyden extrapolation coefficients are not linear in the residual history), where CP2K
//  keeps independent per-channel histories on (ρ,m).  That choice is deliberately left HERE -- swapping the
//  channel basis, or giving m a plain linear leaf while ρ keeps Kerker, is a change to this class alone.
//---------------------------------------------------------------------------------------------------------

//! The polarized ρ̃ Fock density: the two channel mixers' outputs presented as ONE density carrying BOTH
//! faces the framework needs.  Hartree reads the ↑+↓ TOTAL through \c FourierDensity (spin never enters
//! V_H); the spin-native XC engine reads the CHANNELS through \c cSpinResolved_CD -- which is precisely
//! what a single-map ρ̃ mixer could not provide.  A non-owning VIEW: the channel densities belong to the
//! leaf mixers, which outlive every Fock build they feed; \c Seat re-points it after each mix.
class PolarizedMixCD
    : public virtual tChargeDensity<dcmplx>
    , public virtual FourierDensity
    , public virtual cSpinResolved_CD
{
public:
    void Seat(const cChargeDensity* up, const cChargeDensity* dn) { itsUp=up; itsDn=dn; }

    //! cSpinResolved_CD -- the spin-native XC engine's channel access (the whole point of this class).
    virtual const cChargeDensity* GetChannel(const Spin& s) const
    {
        assert(s!=Spin::None && "PolarizedMixCD::GetChannel: ask for a channel, not the total");
        assert(itsUp && itsDn && "PolarizedMixCD: never seated");
        return s==Spin::Up ? itsUp : itsDn;
    }
    // FourierDensity -- the ↑+↓ TOTAL.  Both quantities are LINEAR in ρ̃ (V_H = 4π ρ̃/|G|²), so summing the
    // channels' answers IS the total's answer.
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
    { return MapAdd(Fourier(itsUp).GetFourierDensity(c), Fourier(itsDn).GetFourierDensity(c)); }
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
    { return MapAdd(Fourier(itsUp).GetRepulsion3C(c), Fourier(itsDn).GetRepulsion3C(c)); }
    //! The raw-raster shadow, summed -- empty (pipeline off) unless BOTH channels answer on the same raster.
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
    {
        rvec_t u=Fourier(itsUp).GetRhoOnGrid(c), d=Fourier(itsDn).GetRhoOnGrid(c);
        if (u.size()==0 || u.size()!=d.size()) return rvec_t{};
        u+=d;
        return u;
    }
    // ScalarFunction<double> / tChargeDensity -- the total ρ(r), summed.
    //! Batch evaluation DELEGATES to the channels' batch form.  The inherited default (a loop over the
    //! single-point \c operator() ) would re-walk each channel's ρ̃ map per point and throw away exactly the
    //! flattening + phase-factor hoisting \c FourierMixCD::EvalBatch exists to provide -- and the XC mesh
    //! asks for tens of thousands of points every iteration.
    virtual rvec_t EvalBatch(const rvec3vec_t& r) const
    {
        rvec_t u=itsUp->EvalBatch(r);
        u+=itsDn->EvalBatch(r);
        return u;
    }
    virtual double  operator()(const rvec3_t& r) const { return (*itsUp)(r) + (*itsDn)(r); }
    virtual rvec3_t Gradient  (const rvec3_t& r) const { return itsUp->Gradient(r) + itsDn->Gradient(r); }
    virtual double  GetTotalCharge() const { return itsUp->GetTotalCharge() + itsDn->GetTotalCharge(); }
    //! Read-only VIEW: scaling belongs to the channels the leaf mixers own, never to their presentation.
    virtual void    ReScale(double) { assert(false && "PolarizedMixCD is a read-only view of the mixers' channels"); }
    //! Self-maintaining logical-clock serial.  DERIVED from the channels' serials rather than stamped in
    //! \c Seat on purpose: a mixer frees its old channel density and allocates the new one, which can land at
    //! the SAME address, so pointer identity is not a safe change detector — that is exactly the stale-cache
    //! failure [[project_hamiltonian_dynamic_cache_bug]] cured by moving to serials.  Channel serials come
    //! from the one global counter, so they are unique and never reused.
    virtual size_t Version() const
    {
        const size_t u=itsUp->Version(), d=itsDn->Version();
        if (u!=itsUpV || d!=itsDnV) { itsUpV=u; itsDnV=d; itsVersion=NextDensityVersion(); }
        return itsVersion;
    }
private:
    //! The channel's Fourier face -- a ρ̃ mixer's output always has it (checked at construction, not here).
    static const FourierDensity& Fourier(const cChargeDensity* cd)
    {
        auto* f=dynamic_cast<const FourierDensity*>(cd);
        assert(f && "PolarizedMixCD: a channel density must carry the FourierDensity face");
        return *f;
    }
    const cChargeDensity* itsUp=nullptr;
    const cChargeDensity* itsDn=nullptr;
    mutable size_t itsVersion=0, itsUpV=0, itsDnV=0;   // 0 = the reserved "no density yet" sentinel
};

//! ONE ρ̃ mixer per spin channel, plus the spin-resolved Fock view they feed.  Pure forwarding: the leaves
//! are ordinary \c tDensityMixer<dcmplx>s built by the SAME factory the unpolarized path uses, so Kerker and
//! Pulay (and their successors) are supported here without this class knowing which it holds.
//! WHICH pair of linear combinations the two leaves mix.
enum class ChannelBasis
{
    SpinChannels,    //!< (ρ↑, ρ↓) -- one leaf each.  Reproduces CP2K's Kerker exactly (see the note above).
    TotalAndMoment   //!< (ρ, m=ρ↑−ρ↓) -- lets the two get DIFFERENT leaves.  The physics reason to want it:
                     //!< Kerker's G²/(G²+G₀²) models the Hartree restoring force against long-wavelength
                     //!< CHARGE fluctuations, and the magnetization has no such force, so damping its low-G
                     //!< residual is unmotivated.  Pairing a Kerker leaf on ρ with a PLAIN LINEAR leaf on m
                     //!< (= Kerker at G₀=0) is a 2x2 COUPLED mix in spin space -- provably NOT reachable by
                     //!< any per-(ρ↑,ρ↓) leaf pair, which is why the G₀ sweep could not test it.  Here we
                     //!< deliberately DIVERGE from CP2K, which damps m with the same filter as ρ.
};

class PolarizedDensityMixer : public tDensityMixer<dcmplx>
{
public:
    //! \a up/\a dn are the two leaves.  Under \c SpinChannels they mix ρ↑/ρ↓; under \c TotalAndMoment they
    //! mix ρ/m and MUST both satisfy \c cFourierMixer (checked here) -- m is a difference of maps, so it can
    //! only be driven through the G-space face.  \a fit/\a recip are needed in that mode to read the working
    //! channels' ρ̃ and to rebuild the channel pair the Fock consumes.
    //! (ρ↑,ρ↓) basis: one leaf per spin channel.  Needs nothing but the leaves.
    PolarizedDensityMixer(std::unique_ptr<tDensityMixer<dcmplx>> up,
                          std::unique_ptr<tDensityMixer<dcmplx>> dn)
        : itsUp(std::move(up)), itsDn(std::move(dn)), itsBasis(ChannelBasis::SpinChannels) {}

    //! (ρ,m) basis: \a rho and \a m mix the total and the magnetization, so BOTH must be field mixers
    //! (m is a difference of maps -- see GField).  \a fit/\a recip read the working channels' ρ̃ and give
    //! the rebuilt channel pair its metric.  Separate constructor rather than defaulted arguments: in this
    //! basis they are REQUIRED, and there is no meaningful "no reciprocal lattice".
    PolarizedDensityMixer(std::unique_ptr<tDensityMixer<dcmplx>> rho,
                          std::unique_ptr<tDensityMixer<dcmplx>> m,
                          std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
                          const ReciprocalLattice& recip)
        : itsUp(std::move(rho)), itsDn(std::move(m)), itsBasis(ChannelBasis::TotalAndMoment)
        , itsFit(std::move(fit)), itsRecip(recip)
    {
        itsRhoF=dynamic_cast<tFieldMixer*>(itsUp.get());
        itsMF  =dynamic_cast<tFieldMixer*>(itsDn.get());
        assert(itsRhoF && itsMF && itsFit &&
               "PolarizedDensityMixer(ρ,m): both leaves must be field mixers");
    }

    //! Mix each channel with its own leaf; the SCF gate is the WORSE channel (both must converge -- an AFM
    //! solution whose total has settled while the staggering still moves is not converged).
    double Mix(cd_t& working, const cd_t& old) override
    {
        auto [wu,wd]=Channels(working);
        if (itsBasis==ChannelBasis::SpinChannels)
        {
            auto [ou,od]=Channels(old);
            const double du=itsUp->Mix(wu,ou), dd=itsDn->Mix(wd,od);
            return std::max(du,dd);
        }
        // (ρ,m): read BOTH channels' fresh ρ̃, form the two combinations, and give each its own leaf.
        const FourierDensity& fu=FourierOf(wu.get());
        const FourierDensity& fd=FourierOf(wd.get());
        const ΔG_Map up=fu.GetFourierDensity(*itsFit), dn=fd.GetFourierDensity(*itsFit);
        const rvec_t rup=fu.GetRhoOnGrid(*itsFit),     rdn=fd.GetRhoOnGrid(*itsFit);
        const double dRho=itsRhoF->MixField({MapAdd(up,dn), RawCombine(rup,rdn,+1.0,1.0)});
        const double dM  =itsMF  ->MixField({MapSub(up,dn), RawCombine(rup,rdn,-1.0,1.0)});
        RebuildChannels(wu->GetTotalCharge(), wd->GetTotalCharge());
        return std::max(dRho,dM);
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t& working) const override
    {
        if (itsBasis==ChannelBasis::SpinChannels)
        {
            auto [wu,wd]=Channels(working);
            itsFock.Seat(itsUp->FockDensity(wu), itsDn->FockDensity(wd));
        }
        else if (!itsChUp)                        // iteration 1: FockDensity runs BEFORE the first Mix
        {
            auto [wu,wd]=Channels(working);
            RebuildChannels(wu->GetTotalCharge(), wd->GetTotalCharge());
        }
        return &itsFock;
    }
    double      GetRelax() const override { return itsUp->GetRelax(); }
    const char* Tag     () const override { return itsUp->Tag(); }   // the trace reports the LEAF recipe
    bool   WantsReDamp(const MixSignals& s) const override
    { return itsUp->WantsReDamp(s) || itsDn->WantsReDamp(s); }        // either channel asking is enough
    double ReDampMix(cd_t& working, const cd_t& old) override
    {
        auto [wu,wd]=Channels(working);
        auto [ou,od]=Channels(old);
        return std::max(itsUp->ReDampMix(wu,ou), itsDn->ReDampMix(wd,od));
    }
    void UpdateRelax(const MixSignals& s) override { itsUp->UpdateRelax(s); itsDn->UpdateRelax(s); }

private:
    //! The two spin channels as ALIASING shared_ptrs: they share OWNERSHIP with the parent polarized density
    //! (which therefore cannot die under a leaf) while pointing at the channel that leaf mixes.  A leaf that
    //! mutates its working density (the linear mixer's MixIn) then edits the channel in place, which is
    //! exactly right; the ρ̃ leaves only read it through the FourierDensity face.
    static std::pair<cd_t,cd_t> Channels(const cd_t& cd)
    {
        auto* pol=dynamic_cast<tPolarized_CD<dcmplx>*>(cd.get());
        assert(pol && "PolarizedDensityMixer: the working density must be polarized");
        return { cd_t(cd, pol->GetChargeDensity(Spin::Up)), cd_t(cd, pol->GetChargeDensity(Spin::Down)) };
    }
    static const FourierDensity& FourierOf(const tChargeDensity<dcmplx>* cd)
    {
        auto* f=dynamic_cast<const FourierDensity*>(cd);
        assert(f && "PolarizedDensityMixer: a spin channel must carry the FourierDensity face");
        return *f;
    }
    //! (ρ,m) mode only: rebuild the CHANNEL pair the Fock consumes, ρ̃_σ = (ρ̃ ± m̃)/2, from the two leaves'
    //! mixed combinations, and re-seat the view on them.  Rebuilt (rather than viewed) because a fresh
    //! FourierMixCD is what carries the batch evaluator, the Poisson kernel and a fresh logical-clock serial.
    void RebuildChannels(double qUp, double qDn) const
    {
        const ΔG_Map& r=itsRhoF->Mixed().RhoTilde();
        const ΔG_Map& m=itsMF  ->Mixed().RhoTilde();
        itsChUp=std::make_shared<FourierMixCD>(MapScale(MapAdd(r,m),0.5), *itsRecip, qUp);
        itsChDn=std::make_shared<FourierMixCD>(MapScale(MapSub(r,m),0.5), *itsRecip, qDn);
        const rvec_t rr=itsRhoF->Mixed().GetRhoOnGrid(*itsFit), mm=itsMF->Mixed().GetRhoOnGrid(*itsFit);
        if (rvec_t u=RawCombine(rr,mm,+1.0,0.5); u.size()) itsChUp->SetRawRho(std::move(u));
        if (rvec_t d=RawCombine(rr,mm,-1.0,0.5); d.size()) itsChDn->SetRawRho(std::move(d));
        itsFock.Seat(itsChUp.get(), itsChDn.get());
    }
    std::unique_ptr<tDensityMixer<dcmplx>> itsUp, itsDn;
    ChannelBasis itsBasis;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsFit;   //!< (ρ,m) mode: reads the channels' ρ̃
    std::optional<ReciprocalLattice> itsRecip;             //!< (ρ,m) mode: metric of the rebuilt channels
    tFieldMixer *itsRhoF=nullptr, *itsMF=nullptr;          //!< the leaves' FIELD face ((ρ,m) mode)
    mutable std::shared_ptr<FourierMixCD> itsChUp, itsChDn;//!< the rebuilt channel pair ((ρ,m) mode)
    mutable PolarizedMixCD itsFock;   //!< re-seated by FockDensity (the leaves' outputs change every mix)
};

//! Build ONE ρ̃-space mixer -- Pulay when \a pulayDepth>0, else Kerker -- seeded from \a seed's ρ̃.  The whole
//! recipe in one place, so the polarized path COMPOSES two of these rather than duplicating it.  \a seed is
//! the whole density on the unpolarized path and one spin channel on the polarized one; nothing else differs.
//! \a label names the arm in the one-time banner ("" / "↑" / "↓"), which is emitted HERE because this is
//! where the ρ̃ map and the raw-raster shadow are in hand -- a caller-side banner would rebuild both.
inline std::unique_ptr<tDensityMixer<dcmplx>> MakeGSpaceMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart,
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip,
    GField seed0, double charge, const std::string& label="")
{
    ΔG_Map rho0 = std::move(seed0.tilde);
    rvec_t raw0 = std::move(seed0.raster);       // raw-raster shadow seed (0.5(f2)); empty = late-activate
    std::cerr << "[" << (pulayDepth>0 ? "Pulay" : "Kerker") << "] ENABLED"
              << (label.empty() ? std::string() : " ("+label+")") << ": G0=" << kerkerG0
              << ", ρ̃-mixing on " << charge << " electrons (" << rho0.size() << " G-vectors"
              << (raw0.size() ? ", raw-XC shadow ON" : "") << ")." << std::endl;
    if (pulayDepth>0)
        return std::make_unique<PulayMixer>(relax0, kerkerG0, pulayDepth, pulayStart, std::move(fit), recip,
                                            std::move(rho0), charge, std::move(raw0));
    auto mixed = std::make_shared<FourierMixCD>(std::move(rho0), recip, charge);
    return std::make_unique<KerkerMixer>(relax0, kerkerG0, std::move(fit), std::move(mixed), std::move(raw0));
}

//! Convenience overload seeded from a DENSITY -- it just reads the field off the density's Fourier face.
//! (The field-taking overload above is the primary one: it is what m can be built from.)
inline std::unique_ptr<tDensityMixer<dcmplx>> MakeGSpaceMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart,
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit, const ReciprocalLattice& recip,
    const tChargeDensity<dcmplx>* seed, const std::string& label="")
{
    auto* fd = dynamic_cast<const FourierDensity*>(seed);
    assert(fd && "MakeGSpaceMixer: the seed density must carry the FourierDensity face");
    return MakeGSpaceMixer(relax0, kerkerG0, pulayDepth, pulayStart, fit, recip,
                           GField{fd->GetFourierDensity(*fit), fd->GetRhoOnGrid(*fit)},
                           seed->GetTotalCharge(), label);
}

//! Build the density mixer for a run: Kerker when \a kerkerG0>0 AND the basis/cell/seed are periodic
//! (Band_FT_IBS + UnitCell + FourierDensity), else linear D-mixing.  Mirrors the old KerkerSetup, incl. the
//! LOUD fall-back (Release has no asserts).  \a relax0 = StartingRelaxRo (α=1 => passthrough).
template <class T> std::unique_ptr<tDensityMixer<T>> MakeDensityMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart, const BasisSet::tBasisSet<T>* basis,
    const Structure* structure, const tDM_CD<T>* seed)
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        if (kerkerG0>0.0 || pulayDepth>0)   // both need the periodic ρ̃ machinery
        {
            auto* ftb  = basis ? dynamic_cast<const BasisSet::Band_FT_IBS*>((*basis)[0]) : nullptr;
            auto* cell = dynamic_cast<const UnitCell*>(structure);
            auto* fd   = dynamic_cast<const FourierDensity*>(seed);
            if (!ftb || !cell || !fd)
            {
                std::cerr << "[Mixer] DISABLED: Kerker/Pulay need a periodic Band_FT_IBS basis + UnitCell + "
                          << "FourierDensity -- falling back to linear D-mixing." << std::endl;
                return std::make_unique<LinearMixer<T>>(relax0);
            }
            auto fit = std::shared_ptr<const BasisSet::cFIT_SF_ABS>(ftb->CreateVxcFitBasisSet(cell, qcMesh::MeshParams{}));
            ReciprocalLattice recip(cell->MakeReciprocalCell());
            const char* tag = pulayDepth>0 ? "Pulay" : "Kerker";
            // POLARIZED: one ρ̃ mixer PER SPIN CHANNEL, composed (see PolarizedDensityMixer).  A single-map
            // mixer would hand the Fock a spin-blind total and collapse v_xc to the ζ=0 branch from
            // iteration 1 -- the MnO AFM-II collapse, 2026-08-07.  (QCHEM_SPINBLIND_KERKER=1 takes the
            // single-map path on a polarized density anyway: the A/B valve that re-measures the collapse,
            // and the negative control behind GPW_SCF.PolarizedRunKeepsItsSpin.  Never a production setting.)
            const bool spinBlind = std::getenv("QCHEM_SPINBLIND_KERKER");
            if (auto* pol = spinBlind ? nullptr : dynamic_cast<const tPolarized_CD<T>*>(seed))
            {
                // CHANNEL BASIS.  Default (ρ↑,ρ↓) reproduces CP2K's Kerker exactly.  QCHEM_MIX_RHO_M=1
                // selects (ρ,m): Kerker on ρ, PLAIN LINEAR on m (G0=0 -- the filter is identically 1), the
                // construction the G₀ sweep could NOT test, because a uniform filter on (ρ↑,ρ↓) is
                // algebraically the same operator as that filter on (ρ,m).  Physics: Kerker models the
                // Hartree restoring force against long-wavelength CHARGE fluctuations; m has none.
                if (std::getenv("QCHEM_MIX_RHO_M"))
                {
                    const auto* up=pol->GetChargeDensity(Spin::Up), *dn=pol->GetChargeDensity(Spin::Down);
                    const auto& fu=dynamic_cast<const FourierDensity&>(*up);
                    const auto& fd2=dynamic_cast<const FourierDensity&>(*dn);
                    const ΔG_Map tu=fu.GetFourierDensity(*fit), td=fd2.GetFourierDensity(*fit);
                    const rvec_t ru=fu.GetRhoOnGrid(*fit),      rd=fd2.GetRhoOnGrid(*fit);
                    auto mr=MakeGSpaceMixer(relax0,kerkerG0,pulayDepth,pulayStart,fit,recip,
                                            GField{MapAdd(tu,td), RawCombine(ru,rd,+1.0,1.0)},
                                            up->GetTotalCharge()+dn->GetTotalCharge(), "ρ");
                    auto mm=MakeGSpaceMixer(relax0,/*G0*/0.0,pulayDepth,pulayStart,fit,recip,
                                            GField{MapSub(tu,td), RawCombine(ru,rd,-1.0,1.0)},
                                            0.0, "m: LINEAR, undamped");
                    return std::make_unique<PolarizedDensityMixer>(std::move(mr), std::move(mm), fit, recip);
                }
                auto up=MakeGSpaceMixer(relax0,kerkerG0,pulayDepth,pulayStart,fit,recip,pol->GetChargeDensity(Spin::Up  ),"↑");
                auto dn=MakeGSpaceMixer(relax0,kerkerG0,pulayDepth,pulayStart,fit,recip,pol->GetChargeDensity(Spin::Down),"↓");
                return std::make_unique<PolarizedDensityMixer>(std::move(up), std::move(dn));
            }
            // A spin-resolved density that is NOT a tPolarized_CD cannot hand out mutable channel densities
            // for the leaves to mix -- never silently unpolarized, so say so and take the physical path.
            if (!spinBlind && dynamic_cast<const tSpinResolved_CD<T>*>(seed))
            {
                std::cerr << "[Mixer] " << tag << " DISABLED: this spin-resolved density exposes no mutable "
                          << "channels to mix per spin -- falling back to linear D-mixing (which keeps both "
                          << "channels) rather than collapsing v_xc to the unpolarized branch." << std::endl;
                return std::make_unique<LinearMixer<T>>(relax0);
            }
            return MakeGSpaceMixer(relax0, kerkerG0, pulayDepth, pulayStart, fit, recip, seed);
        }
    }
    return std::make_unique<LinearMixer<T>>(relax0);
}

} //namespace
