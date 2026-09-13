// File: ChargeDensity/Internal/PolarizedDensityMixer.C  The polarized ρ̃ mixer: one G-space leaf per channel
// ((ρ↑,ρ↓) or (ρ,m)), one JOINT history, and the spin-resolved Fock presentation they feed.
module;
#include <memory>
#include <iostream>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <vector>
export module qchem.ChargeDensity.Internal.PolarizedDensityMixer;
export import qchem.ChargeDensity.DensityMixer;
export import qchem.ChargeDensity.Internal.FieldMixer;
import qchem.ReciprocalLattice;
import qchem.Blaze;
import qchem.Types;

export namespace qchem::ChargeDensity
{

//---------------------------------------------------------------------------------------------------------
//  THE POLARIZED ρ̃ MIXER (2026-08-07; doc/SymmetryUpgradePlan.md §7 step 7)
//
//  A ρ̃-space mixer carries ONE FourierMixCD -- the ↑+↓ total, no spin channels -- and drives every Fock
//  from it.  DensitySampler::RhoPol then finds no spin face and takes its ρ↑=ρ↓=ρ/2 branch, so v_xc^↑≡v_xc^↓
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
//  mixers (Pulay/Broyden extrapolation coefficients are not linear in the residual history).  That choice is
//  deliberately left HERE -- swapping the channel basis, or giving m a plain linear leaf while ρ keeps
//  Kerker, is a change to this class alone.
//
//  WHAT THE HISTORY MAY *NOT* DO (2026-08-10, the MnO ejections; §7 step 7).  Splitting the FILTER per
//  channel is free; splitting the HISTORY is not.  This class ran ONE Pulay per channel, i.e. two
//  independent B-solves giving c↑ ≠ c↓, so the extrapolated state (Σcᵢ↑ρ↑ᵢ, Σcᵢ↓ρ↓ᵢ) never occurred on the
//  trajectory: each channel still conserved its charge, but the MOMENT came out an arbitrary synthesised
//  combination of history moments -- MnO run 27's ejections (E −61.25 → −49, m_stag 0.4 → 0.1, recovering,
//  repeating).  So this class is now a channel-basis PRECONDITIONER feeding ONE joint extrapolator:
//      residual → per-channel filter (may differ per channel) → joint extrapolation (one B, one c, all channels)
//  which is also the VASP/QE/CP2K architecture -- and note CP2K's BETA 1.5 IS its Kerker preconditioner
//  sitting in front of ONE Broyden history.  The two concepts COMPOSE; we had them fused.
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
    , public virtual cDM_SourceSink   //!< the driver seats the polarized D here; the view splits it per channel
{
public:
    void Seat(const cChargeDensity* up, const cChargeDensity* dn) { itsUp=up; itsDn=dn; }

    //! cDM_SourceSink -- split the polarized D into its channels and seat each on the channel it belongs to.
    //! ALIASING shared_ptrs: each channel pointer keeps the PARENT density alive while pointing at the child,
    //! the ownership the channel accessors (raw, non-owning) cannot express on their own -- and a retained
    //! XC source genuinely outlives the call, unlike a mix.  A non-polarized D (the iteration-0 seed) seats
    //! nothing.  Seated on whatever channels the view CURRENTLY presents -- which in (ρ,m) mode are the
    //! rebuilt ones, closing a hole the mixer-side deposit had: it reached the leaves, never the rebuilt pair.
    virtual void SetDMSource(std::shared_ptr<const cDM_CD> dm) const override
    {
        auto* pol = dynamic_cast<const cPolarized_CD*>(dm.get());
        if (!pol) return;
        Sink(itsUp).SetDMSource(std::shared_ptr<const cDM_CD>(dm, pol->GetChargeDensity(Spin::Up  )));
        Sink(itsDn).SetDMSource(std::shared_ptr<const cDM_CD>(dm, pol->GetChargeDensity(Spin::Down)));
    }

    //! cSpinResolved_CD -- the spin-native XC engine's channel access (the whole point of this class).
    virtual const cChargeDensity* GetChannel(const Spin& s) const override
    {
        assert(s!=Spin::None && "PolarizedMixCD::GetChannel: ask for a channel, not the total");
        assert(itsUp && itsDn && "PolarizedMixCD: never seated");
        return s==Spin::Up ? itsUp : itsDn;
    }
    // FourierDensity -- the ↑+↓ TOTAL.  Both quantities are LINEAR in ρ̃ (V_H = 4π ρ̃/|G|²), so summing the
    // channels' answers IS the total's answer.
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const override
    { return MapAdd(Fourier(itsUp).GetFourierDensity(c), Fourier(itsDn).GetFourierDensity(c)); }
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const override
    { return MapAdd(Fourier(itsUp).GetRepulsion3C(c), Fourier(itsDn).GetRepulsion3C(c)); }
    //! The raw-raster shadow, summed -- empty (pipeline off) unless BOTH channels answer on the same raster.
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const override
    {
        rvec_t u=Fourier(itsUp).GetRhoOnGrid(c), d=Fourier(itsDn).GetRhoOnGrid(c);
        if (u.size()==0 || u.size()!=d.size()) return rvec_t{};
        u+=d;
        return u;
    }
    // ScalarFunction<double> / tChargeDensity -- the total ρ(r), summed.
    //! Batch evaluation DELEGATES to the channels' batch form.  The inherited default (a loop over the
    //! single-point \c operator() ) would re-walk each channel's ρ̃ map per point and throw away exactly the
    //! flattening + phase-factor hoisting \c FourierMixCD's batch override exists to provide -- and the XC
    //! mesh asks for tens of thousands of points every iteration.
    //! (Spelled as the inherited \c ScalarFunction batch operator: the parallel \c EvalBatch fork this used
    //! to call was deleted by R1.5, precisely so a neutral-face caller cannot miss the fast path.)
    virtual rvec_t operator()(const rvec3vec_t& r) const override
    {
        rvec_t u=(*itsUp)(r);
        u+=(*itsDn)(r);
        return u;
    }
    virtual double  operator()(const rvec3_t& r) const override { return (*itsUp)(r) + (*itsDn)(r); }
    virtual rvec3_t Gradient  (const rvec3_t& r) const override { return itsUp->Gradient(r) + itsDn->Gradient(r); }
    virtual double  GetTotalCharge() const override { return itsUp->GetTotalCharge() + itsDn->GetTotalCharge(); }
    //! Read-only VIEW: scaling belongs to the channels the leaf mixers own, never to their presentation.
    virtual void    ReScale(double) override { assert(false && "PolarizedMixCD is a read-only view of the mixers' channels"); }
    //! Self-maintaining logical-clock serial.  DERIVED from the channels' serials rather than stamped in
    //! \c Seat on purpose: a mixer frees its old channel density and allocates the new one, which can land at
    //! the SAME address, so pointer identity is not a safe change detector — that is exactly the stale-cache
    //! failure [[project_hamiltonian_dynamic_cache_bug]] cured by moving to serials.  Channel serials come
    //! from the one global counter, so they are unique and never reused.
    virtual size_t Version() const override
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
    //! The channel's sink face -- a ρ̃ mixer's output always has it (a FourierMixCD).
    static const cDM_SourceSink& Sink(const cChargeDensity* cd)
    {
        auto* k=dynamic_cast<const cDM_SourceSink*>(cd);
        assert(k && "PolarizedMixCD: a channel density must carry the DM-source sink face");
        return *k;
    }
    const cChargeDensity* itsUp=nullptr;
    const cChargeDensity* itsDn=nullptr;
    mutable size_t itsVersion=0, itsUpV=0, itsDnV=0;   // 0 = the reserved "no density yet" sentinel
};

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

//! \brief The CHANNEL-BASIS PRECONDITIONER for a polarized ρ̃ run: two leaves that FILTER a pair of channels
//! -- (ρ↑,ρ↓) or (ρ,m) -- plus the spin-resolved Fock view they feed.  When the leaves carry HISTORY it does
//! NOT let them extrapolate independently: it stages both and runs ONE joint solve (see \c MixJointly).
//!
//! The leaves are ordinary mixers built by the SAME factory the unpolarized path uses, so Kerker and Pulay
//! (and their successors) are supported without this class knowing which it holds -- but that ignorance now
//! stops exactly where it must: it ASKS each leaf whether it carries memory (\c tFieldMixer::History), because
//! that is precisely the property deciding whether splitting a step per channel is legitimate.
class PolarizedDensityMixer : public tDensityMixer<dcmplx>
{
public:
    //! \a a/\a b are the two leaves and must BOTH be field mixers: under \c SpinChannels they mix ρ↑/ρ↓,
    //! under \c TotalAndMoment ρ and m=ρ↑−ρ↓ (m is a difference of maps, so it can only be driven through the
    //! G-space face -- see \c GField).  \a fit reads the working channels' ρ̃; \a recip gives the rebuilt
    //! channel pair its metric.  Both are REQUIRED in either basis: the channel fields are formed HERE now,
    //! so that a history-carrying pair can be staged before either is committed.
    PolarizedDensityMixer(std::unique_ptr<tDensityMixer<dcmplx>> a,
                          std::unique_ptr<tDensityMixer<dcmplx>> b,
                          std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
                          const ReciprocalLattice& recip,
                          ChannelBasis basis)
        : itsUp(std::move(a)), itsDn(std::move(b)), itsBasis(basis)
        , itsFit(std::move(fit)), itsRecip(recip)
    {
        itsAF=dynamic_cast<tFieldMixer*>(itsUp.get());
        itsBF=dynamic_cast<tFieldMixer*>(itsDn.get());
        assert(itsAF && itsBF && itsFit && "PolarizedDensityMixer: both leaves must be field mixers");
        itsAH=itsAF->History(); itsBH=itsBF->History();
        if (itsAH && itsBH)
            std::cerr << "[Pulay] JOINT HISTORY: one B summed over both channels, one coefficient vector "
                      << "(spin is just another irrep)." << std::endl;
    }

    //! Read both channels' fresh ρ̃, form the pair of channel fields, and mix them.  The SCF gate is the WORSE
    //! channel (both must converge -- an AFM solution whose total has settled while the staggering still moves
    //! is not converged).
    double Mix(cd_t& working, const cd_t& /*old*/) override
    {
        auto [wu,wd]=Channels(working);
        const FourierDensity& fu=FourierOf(wu);
        const FourierDensity& fd=FourierOf(wd);
        const ΔG_Map up=fu.GetFourierDensity(*itsFit), dn=fd.GetFourierDensity(*itsFit);
        const rvec_t rup=fu.GetRhoOnGrid(*itsFit),     rdn=fd.GetRhoOnGrid(*itsFit);
        const bool spin=(itsBasis==ChannelBasis::SpinChannels);
        const GField fa = spin ? GField{up, rup} : GField{MapAdd(up,dn), RawCombine(rup,rdn,+1.0,1.0)};
        const GField fb = spin ? GField{dn, rdn} : GField{MapSub(up,dn), RawCombine(rup,rdn,-1.0,1.0)};
        // Every HISTORY-carrying channel goes into ONE extrapolation; a memoryless leaf just filters its own
        // channel, which for a channel-diagonal filter is the same operator either way.  (All four leaf
        // combinations are legitimate -- history on ρ with a plain filter on m is a real recipe -- and none
        // of them splits a history, which is the only thing forbidden.)
        std::vector<tFieldExtrapolator*> hist;
        std::vector<GField>              histFields;
        double d=0.0;
        if (itsAH) { hist.push_back(itsAH); histFields.push_back(fa); } else d=std::max(d,itsAF->MixField(fa));
        if (itsBH) { hist.push_back(itsBH); histFields.push_back(fb); } else d=std::max(d,itsBF->MixField(fb));
        if (!hist.empty()) d=std::max(d,MixJointly(hist,histFields));
        if (!spin) RebuildChannels(wu->GetTotalCharge(), wd->GetTotalCharge());
        return d;
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t& working) const override
    {
        if (itsBasis==ChannelBasis::SpinChannels)
            itsFock.Seat(&itsAF->Mixed(), &itsBF->Mixed());   // re-seat: each mix allocates a fresh mixed ρ̃
        else if (!itsChUp)                        // iteration 1: FockDensity runs BEFORE the first Mix
        {
            auto [wu,wd]=Channels(working);
            RebuildChannels(wu->GetTotalCharge(), wd->GetTotalCharge());
        }
        return &itsFock;
    }
    double      GetRelax() const override { return itsUp->GetRelax(); }
    const char* Tag     () const override { return itsUp->Tag(); }   // the trace reports the LEAF recipe
    // (No adaptive hooks: the leaves are G-space mixers, none of which adapts its step -- V1.18.)

private:
    //! The two spin channels.  Plain pointers now that the mixer's subject is a reference: the parent
    //! density is owned by the caller and outlives every call, so the aliasing shared_ptrs this used to
    //! build (purely to keep the parent alive under a leaf) are gone with the shared_ptr itself.
    //! A leaf that mutates its working density edits the channel in place, which is exactly right.
    static std::pair<cd_t*,cd_t*> Channels(cd_t& cd)
    {
        auto* pol=dynamic_cast<tPolarized_CD<dcmplx>*>(&cd);
        assert(pol && "PolarizedDensityMixer: the working density must be polarized");
        return { pol->GetChargeDensity(Spin::Up), pol->GetChargeDensity(Spin::Down) };
    }
    //! const overload -- the \a old density is only ever READ (constness carried through, no const_cast).
    static std::pair<const cd_t*,const cd_t*> Channels(const cd_t& cd)
    {
        auto* pol=dynamic_cast<const tPolarized_CD<dcmplx>*>(&cd);
        assert(pol && "PolarizedDensityMixer: the working density must be polarized");
        return { pol->GetChargeDensity(Spin::Up), pol->GetChargeDensity(Spin::Down) };
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
        const ΔG_Map& r=itsAF->Mixed().RhoTilde();
        const ΔG_Map& m=itsBF->Mixed().RhoTilde();
        itsChUp=std::make_shared<FourierMixCD>(MapScale(MapAdd(r,m),0.5), itsRecip, qUp);
        itsChDn=std::make_shared<FourierMixCD>(MapScale(MapSub(r,m),0.5), itsRecip, qDn);
        const rvec_t rr=itsAF->Mixed().GetRhoOnGrid(*itsFit), mm=itsBF->Mixed().GetRhoOnGrid(*itsFit);
        if (rvec_t u=RawCombine(rr,mm,+1.0,0.5); u.size()) itsChUp->SetRawRho(std::move(u));
        if (rvec_t d=RawCombine(rr,mm,-1.0,0.5); d.size()) itsChDn->SetRawRho(std::move(d));
        itsFock.Seat(itsChUp.get(), itsChDn.get());
    }
    std::unique_ptr<tDensityMixer<dcmplx>> itsUp, itsDn;
    ChannelBasis itsBasis;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsFit;   //!< reads the working channels' ρ̃
    ReciprocalLattice itsRecip;                            //!< metric of the rebuilt channels ((ρ,m) mode)
    tFieldMixer *itsAF=nullptr, *itsBF=nullptr;            //!< the leaves' FIELD face
    tFieldExtrapolator *itsAH=nullptr, *itsBH=nullptr;     //!< ...and their HISTORY face (null = memoryless)
    mutable std::shared_ptr<FourierMixCD> itsChUp, itsChDn;//!< the rebuilt channel pair ((ρ,m) mode)
    mutable PolarizedMixCD itsFock;   //!< re-seated by FockDensity (the leaves' outputs change every mix)
};

} //namespace
