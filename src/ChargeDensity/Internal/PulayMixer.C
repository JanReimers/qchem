// File: ChargeDensity/Internal/PulayMixer.C  Density-DIIS (Pulay) ρ̃ mixing, Kerker-preconditioned.
module;
#include <memory>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <deque>
#include <vector>
export module qchem.ChargeDensity.Internal.PulayMixer;
export import qchem.ChargeDensity.DensityMixer;
export import qchem.ChargeDensity.Internal.FieldMixer;
import qchem.BasisSet.G_FieldEvaluator;            // G_SpectralFilter (the raster shadow's evaluator)
import qchem.ReciprocalLattice;
import qchem.Blaze;
import qchem.Types;

export namespace qchem::ChargeDensity
{

//! PULAY (density-DIIS) ρ̃-mixing, Kerker-preconditioned (periodic / dcmplx).  Keeps a history of the fed
//! densities ρ̃_in and the freshly collocated ρ̃_out; each step solves the DIIS bordered system (the shared
//! qchem.Math.DIIS engine) over the residuals ρ̃_out−ρ̃_in for the optimal coefficients c (Σc=1), forms the
//! extrapolated ρ̃_in*=Σcᵢρ̃_inᵢ and ρ̃_out*=Σcᵢρ̃_outᵢ, and applies the Kerker step to THOSE (KerkerStep:
//! ρ̃_in* + α·G²/(G²+G0²)·(ρ̃_out*−ρ̃_in*)).  First iteration (history<2) falls back to plain
//! Kerker.  doc/SCFStrategyPlan.md §4 (the density-face use of the shared extrapolator).
//!
//! On a POLARIZED run one of these mixes each channel, but they SHARE one extrapolation: the step is split
//! into \c StageResidual (grow the history, hand out this channel's block of B) and \c ApplyJoint (commit
//! with the coefficients solved from the summed B) -- see \c tFieldExtrapolator.  The Kerker step itself
//! stays per-channel, because a FILTER may legitimately differ per channel (that is what the (ρ,m) basis
//! buys); only the HISTORY must be shared.
class PulayMixer : public tDensityMixer<dcmplx>, public virtual tFieldExtrapolator
{
public:
    //! \a seed0 is the running field's starting value (read off the seed density by the factory); presented
    //! at once so iteration 0's Fock has a density to be driven from.
    PulayMixer(double relax, double G0, int depth, int start, std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit,
               ReciprocalLattice recip, GField seed0, double charge)
        : itsRelax(relax), itsKerkerG0(G0), itsDepth(depth), itsStart(start), itsKerkerFit(std::move(fit))
        , itsRecip(recip), itsCharge(charge), itsIn(std::move(seed0))
        , itsMixedRho(Present(itsIn, itsRecip, itsCharge))
    {}

    //! The density-face entry: extract ρ̃_out + the raster shadow, then do the G-space arithmetic below.
    double Mix(cd_t& working, const cd_t&) override
    {
        auto* fd = dynamic_cast<const FourierDensity*>(&working);
        assert(fd);
        return MixField({fd->GetFourierDensity(*itsKerkerFit), fd->GetRhoOnGrid(*itsKerkerFit)});
    }
    const GField&       Field() const override { return itsIn; }
    const FourierMixCD& Mixed() const override { assert(itsMixedRho); return *itsMixedRho; }

    //! The single-channel step: stage, solve MY OWN B, apply.  Exactly \c MixJointly over one channel -- and
    //! written as such, so the unpolarized path cannot drift from the polarized one.
    double MixField(const GField& field) override
    {
        return MixJointly({this},{field});
    }

    //! Phase 1 (see \c tFieldExtrapolator): grow the history and hand out this channel's block of B.  Nothing
    //! is committed here -- the pending (in,out) pair is held for \c ApplyJoint.
    StagedResidual StageResidual(const GField& field) override
    {
        assert(!itsStaged && "PulayMixer: StageResidual twice with no ApplyJoint between");
        const ΔG_Map& out    = field.tilde;
        const rvec_t& rawOut = field.raster;
        ΔG_Map in  = itsIn.tilde;                           // ρ̃_in : the field fed to this iteration's Fock
        ΔG_Map res = out-in;                                // residual = ρ̃_out − ρ̃_in
        StagedResidual sr;
        sr.resid = MaxAbs(res);
        // RAW-raster shadow inputs (0.5(f2)); late-activates like KerkerMixer, drops out if answers stop.
        const bool raw = rawOut.size() && RasterEvaluator();
        if (raw && itsIn.raster.size()!=rawOut.size()) itsIn.raster=rawOut;            // bootstrap/late-activate
        if (!raw) { itsIn.raster=rvec_t{}; itsRawIns.clear(); itsRawOuts.clear(); }
        itsPending = Pending{in, out, itsIn.raster, rawOut, raw};   // rawIn = the shadow of `in`, before the update
        itsStaged  = true;

        // PRIME with plain Kerker until we are near the fixed point (history-based mixing is unstable far
        // out).  No history is accumulated during priming, so Pulay starts with clean, linear-regime residuals.
        if (++itsCount<=itsStart) return sr;                // no history ⇒ no B ⇒ ApplyJoint takes the plain step

        itsIns.push_back(in); itsOuts.push_back(out); itsResiduals.push_back(std::move(res));  // grow the history
        if (raw) { itsRawIns.push_back(itsPending.rawIn); itsRawOuts.push_back(rawOut); }
        while ((int)itsResiduals.size()>itsDepth)                                    // prune to `depth`
            { itsIns.pop_front(); itsOuts.pop_front(); itsResiduals.pop_front(); }
        while (itsRawIns.size()>itsResiduals.size()) { itsRawIns.pop_front(); itsRawOuts.pop_front(); }

        const size_t n=itsResiduals.size();
        if (n>=2)                                           // n<2: not enough history yet → plain Kerker
        {
            sr.B=blazem::zero<double>(n);                   // Bᵢⱼ = ⟨resᵢ,resⱼ⟩ (symmetric)
            for (size_t i=0;i<n;++i)
                for (size_t j=i;j<n;++j) sr.B(i,j)=ResidualInnerRe(itsResiduals[i],itsResiduals[j]);
        }
        return sr;
    }

    //! Phase 2: extrapolate with the SHARED coefficients \a c (empty = plain step) and take the Kerker step
    //! on the extrapolated pair.
    void ApplyJoint(const rvec_t& c) override
    {
        assert(itsStaged && "PulayMixer::ApplyJoint with nothing staged");
        itsStaged=false;
        const Pending& p=itsPending;
        ΔG_Map inStar, outStar;
        rvec_t rawInStar=p.rawIn, rawOutStar=p.rawOut;     // raw stars default to the plain pair
        if (c.size()==0) { inStar=p.in; outStar=p.out; }   // priming / no history → plain Kerker
        else
        {
            assert(c.size()==itsResiduals.size() && "PulayMixer::ApplyJoint: c does not span my history");
            inStar =LinearCombination(itsIns ,c);          // ρ̃_in*  = Σ cᵢ ρ̃_inᵢ
            outStar=LinearCombination(itsOuts,c);          // ρ̃_out* = Σ cᵢ ρ̃_outᵢ
            // The raw shadow takes the SAME extrapolation coefficients -- but only once its history spans the
            // whole window (a late-activated shadow falls back to the plain pair until the deques align).
            if (p.raw && itsRawIns.size()==c.size())
            {
                rawInStar =c[0]*itsRawIns[0];  for (size_t i=1;i<c.size();++i) rawInStar +=c[i]*itsRawIns[i];
                rawOutStar=c[0]*itsRawOuts[0]; for (size_t i=1;i<c.size();++i) rawOutStar+=c[i]*itsRawOuts[i];
            }
        }
        // THE MIXER OWNS ITS FIELD (V1.18): the Kerker step on the DIIS-extrapolated pair updates itsIn, and
        // the density the Fock sees is a presentation of the result.
        KerkerStepResult r = KerkerStep(inStar, outStar, itsRelax, itsKerkerG0, itsRecip);
        itsIn.tilde = std::move(r.mix);
        if (p.raw) itsIn.raster = RasterKerker(*RasterEvaluator(), rawInStar, rawOutStar, itsRelax, itsKerkerG0);
        itsMixedRho = Present(itsIn, itsRecip, itsCharge, r.corr, r.alphaEff);
    }
    const tChargeDensity<dcmplx>* FockDensity(const cd_t&) const override { return itsMixedRho.get(); }
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "Pul"; }
private:
    //! The raster-shadow evaluator, or nullptr when the fit basis cannot raster (⇒ the raw pipeline is off).
    const BasisSet::G_SpectralFilter* RasterEvaluator() const
    { return dynamic_cast<const BasisSet::G_SpectralFilter*>(itsKerkerFit.get()); }

    //! The (in,out) pair staged by \c StageResidual and consumed by \c ApplyJoint -- the state the two-phase
    //! split needs and the old single-phase \c MixField kept in locals.
    struct Pending { ΔG_Map in, out; rvec_t rawIn, rawOut; bool raw=false; };

    double itsRelax, itsKerkerG0; int itsDepth, itsStart, itsCount=0;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsKerkerFit;
    ReciprocalLattice itsRecip;
    double itsCharge;
    GField itsIn;                                          //!< THE running mixed field (+ its raster shadow; empty = off)
    std::shared_ptr<FourierMixCD> itsMixedRho;             //!< its presentation, rebuilt each step
    std::deque<ΔG_Map> itsIns, itsOuts, itsResiduals;      // the Pulay history (aligned index-wise)
    std::deque<rvec_t> itsRawIns, itsRawOuts;              //!< its history (aligned with itsIns while active)
    Pending itsPending;                                    //!< staged, not yet committed
    bool    itsStaged=false;
};

} //namespace
