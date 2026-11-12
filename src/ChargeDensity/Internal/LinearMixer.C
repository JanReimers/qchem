// File: ChargeDensity/Internal/LinearMixer.C  The density-MATRIX mixers: linear (legacy adaptive α) and D-space Pulay.
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <cassert>
export module qchem.ChargeDensity.Internal.LinearMixer;
export import qchem.ChargeDensity.DensityMixer;
import qchem.Math.DIIS;   // Bordered / Coefficients -- the shared extrapolator (the ρ̃ Pulay and the Fock DIIS use it too)
import qchem.Blaze;

export namespace qchem::ChargeDensity
{

//! \brief The mixed D_in a density-MATRIX mixer owns (2026-11-09).  The loop's working density is the idempotent D_out of the
//! diagonalisation; the mixed D_in that drives the next Fock is the mixer's own copy, as a field mixer's running field always was.
//! Mixing in place put the mixed, NON-idempotent D in the loop's hands, and everything read it there: the energy pass (first
//! order in the mixing error: the ±1e-8..1e-7 Ha floor of OpenWork §3, measured on Deck.ASavedState and the PW resume), the saved
//! state, the converged density.  \c InputFor guards the hand-back: my D_in drives the Fock only while the working density is the
//! D_out I last mixed (the same \c Version); anything that moved the working density without me (iteration 1, a GDM stretch, a
//! restart) gets itself back, which is what an in-place mixer gave.
template <class T> class OwnedDIn
{
protected:
    typedef tMixableDensity<T> cd_t;
    const cd_t& InputFor(const cd_t& working) const
    { return itsMixed && working.Version()==itsOutVersion ? *itsMixed : working; }
    void Install(std::unique_ptr<cd_t> mixed, const cd_t& out) { itsMixed=std::move(mixed); itsOutVersion=out.Version(); }
    std::unique_ptr<cd_t> itsMixed;    //!< the next Fock's D_in (null: the working D itself)
    size_t                itsOutVersion=0;   //!< the D_out it was mixed from
};

//! Linear density-matrix mixing with the legacy adaptive-α heuristics.  α=1 (the default) = passthrough.
template <class T> class LinearMixer : private OwnedDIn<T>, public tDensityMixer<T>, public virtual tAdaptiveMixer<T>
{
    using OwnedDIn<T>::InputFor; using OwnedDIn<T>::Install;
public:
    typedef typename tDensityMixer<T>::cd_t cd_t;
    //! \a relaxMax: the CEILING the adaptive step may grow to (SCFParams::MaxRelaxRo; a restart recipe sets 0.3).
    //! \a relaxMin: the FLOOR a retraction may shrink it to (SCFParams::MinRelaxRo; 0 = none).  Both are deck keys, so a bad
    //! pair is refused with a throw, not an assert (dead under NDEBUG).
    LinearMixer(double relax0, double relaxMax, double relaxMin=0.0)
        : itsRelax(std::max(std::min(relax0, relaxMax), relaxMin)), itsRelMax(relaxMax), itsRelMaxStated(relaxMax), itsRelMin(relaxMin)
    {
        assert(relaxMax>0.0 && relaxMax<=1.0);
        if (!(relaxMin>=0.0 && relaxMin<=relaxMax))
            throw std::invalid_argument("LinearMixer: the step floor MinRelaxRo ("+std::to_string(relaxMin)+") must lie in [0, "
                                        "MaxRelaxRo = "+std::to_string(relaxMax)+"]");
    }

    //! \a working is this step's D_out and is NOT touched (2026-11-09): the mixed D_in of the next Fock is MINE, so the loop's working
    //! density stays the idempotent D_out every reader wants -- the energy pass, the saved state, the converged density.
    double Mix(const cd_t& working, const cd_t& old) override
    {
        const cd_t& in=InputFor(old);                                          // the D_in that drove THIS Fock
        double dcd = working.GetChangeFrom(in)/working.GetTotalCharge();       // relative change
        if (dcd<1e-5) itsRelMax=std::max(std::min(itsRelMax, 0.5), itsRelMin);   // the floor outranks the late cap
        itsExtrapolated=false;
        itsIn=in.MakeCopy();                                                   // Adapt's retraction rebuilds from it
        itsStepRelax=itsRelax;
        if (itsRelax==1.0) { this->itsMixed.reset(); return dcd; }                  // passthrough: the next Fock IS D_out
        auto mixed=working.MakeCopy();
        mixed->MixIn(in, 1.0-itsRelax);                                        // (1−relax)ρ_in + relax ρ_out
        Install(std::move(mixed), working);
        return dcd;
    }
    //! The accelerator extrapolated this step: a step at my CURRENT α with no late cap, and no retraction can follow (\c Adapt only
    //! grows α toward the STATED ceiling after such a step).  \c tAdaptiveMixer::ExtrapolatedStep has the measurement.
    double ExtrapolatedStep(const cd_t& working, const cd_t& old) override
    {
        const cd_t& in=InputFor(old);                                          // the D_in that drove THIS Fock
        const double dcd = working.GetChangeFrom(in)/working.GetTotalCharge();
        itsIn.reset();                                                         // nothing to retract
        itsExtrapolated=true;
        itsStepRelax=itsRelax;
        if (itsRelax==1.0) { this->itsMixed.reset(); return dcd; }             // the next Fock IS D_out
        auto mixed=working.MakeCopy();
        mixed->MixIn(in, 1.0-itsRelax);
        Install(std::move(mixed), working);
        return dcd;
    }
    const tChargeDensity<T>* FockDensity(const cd_t& working) const override { return &InputFor(working); }
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "Lin"; }
    bool HasGFilter() const override { return false; }   // mixes D itself
    bool PassesOutputThrough() const override { return itsRelax==1.0; }   // Mix's own passthrough test

    //! The [F,D]-keyed policy, in one place: if this step made [F,D] WORSE, retract it to a quarter step and shrink α by 0.8;
    //! if it improved, grow α by the GOLDEN RATIO -- and approach the ceiling without ever landing on it.
    //!
    //! ⚠ WHY IRRATIONAL (user, 2026-10-29: "the ratio of shrink and grow factors should be an irrational number").  With
    //! commensurate factors and a hard clamp the step RETURNS to exactly the same α -- 1.0 x0.8 x1.5 -> clamp 1.0 -- and the
    //! density map follows it into an EXACT cycle: the Si frozen-χ ±α restarts sat on a period-3 limit cycle for 200
    //! iterations (Null accelerator, the moment its [F,D] became real).  ln(phi)/ln(1/0.8) is irrational, and a step that
    //! would cross the ceiling closes only 1/phi^2 of the remaining gap, so α never repeats.  (Measured alternatives: a
    //! ceiling at the last retracted α broke 43 atom/molecule tests; sqrt(2) growth with a 1/phi approach left two too slow.)
    //!
    //! THE RETRACTION is rebuilt from the step's own pair: the quarter step (1−α')D_in + α'D_out, α' = α/4, from the D_in I kept
    //! and \a out (the loop's working D_out, which the retraction leaves alone, as Mix does).
    bool Adapt(const MixSignals& s, const cd_t& out) override
    {
        bool remixed=false;
        constexpr double phi=1.6180339887498949;   // (1+sqrt 5)/2
        if (itsExtrapolated)   // after a step DIIS extrapolated: GROWTH ONLY, toward the stated ceiling -- no retraction, no shrink
        {
            if (s.FD<s.FDold) itsRelax = phi*itsRelax<itsRelMaxStated ? phi*itsRelax : itsRelMaxStated-(itsRelMaxStated-itsRelax)/(phi*phi);
            if (itsRelax>itsRelMaxStated) itsRelax=itsRelMaxStated;
            return false;
        }
        const double dFD = s.FD - s.FDold;
        if (s.FD>s.FDold && std::fabs(dFD)>1e-9 && itsRelax>0.0)
        {
            if (itsIn)   // a step to retract (none before the first Mix)
            {
                auto mixed=itsIn->MakeCopy();
                mixed->MixIn(out, 0.25*itsStepRelax);
                Install(std::move(mixed), out);
                remixed=true;
            }
            itsRelax=std::max(0.8*itsRelax, itsRelMin);
        }
        if (s.FD<s.FDold) itsRelax = phi*itsRelax<itsRelMax ? phi*itsRelax : itsRelMax-(itsRelMax-itsRelax)/(phi*phi);
        if (itsRelax>itsRelMax) itsRelax=itsRelMax;
        return remixed;
    }
private:
    double itsRelax;
    double itsRelMax;   //!< the ceiling (SCFParams::MaxRelaxRo; lowered to 0.5 once the density is nearly converged)
    double itsRelMaxStated;   //!< the ceiling as STATED (the late cap lowers itsRelMax, never this): growth's bound after an extrapolated step
    bool   itsExtrapolated=false;   //!< the last step was \c ExtrapolatedStep: \c Adapt only grows
    double itsRelMin;   //!< the floor (SCFParams::MinRelaxRo): a retraction never shrinks the step below it
    double itsStepRelax=1.0;                 //!< the α the last Mix used (a retraction quarters THAT step)
    std::unique_ptr<tMixableDensity<T>> itsIn;   //!< the D_in of the last Mix (the retraction's base)
};

//! \brief PULAY (density-DIIS) ON THE DENSITY MATRIX -- "Pulay without Kerker" (OpenWork §3, 2026-10-26).
//!
//! The history is D itself: copies of each step's D_in and D_out.  The step solves the DIIS bordered system over the
//! residuals R_i = D_out,i - D_in,i (inner product \c tMixableDensity::Dot) for c (Σc = 1), and the next input is
//! \f$D_{in}^{next}=(1-\alpha)\sum_i c_iD_{in,i}+\alpha\sum_i c_iD_{out,i}\f$, formed on the mixer's OWN copy (\c OwnedDIn; the working D_out is untouched) with
//! \c ReScale + \c Axpy.  So the Fock -- and V_xc -- see that mixed density MATRIX through the DM route: no fit, no ρ̃, no
//! raster shadow, no cusp deficit (pin 28 holds by construction; the loop's RequirePin28 would throw otherwise).  ρ[D] is
//! linear in D, so this IS Pulay on ρ[D] -- exactly, with no G-ball truncation.  Any T, any system: it asks nothing of the
//! geometry.  On a polarized run the composite holds both spins, so ONE history and one B summed over them come free.
//! Priming: \a start plain linear steps at \a relax before the history engages (history-based mixing is unstable far out).
//! Memory: 2 copies of D per history slot (the residual inner products are formed bilinearly from in/out).
template <class T> class PulayDMixer : private OwnedDIn<T>, public tDensityMixer<T>
{
    using OwnedDIn<T>::InputFor; using OwnedDIn<T>::Install;
public:
    typedef typename tDensityMixer<T>::cd_t cd_t;
    PulayDMixer(double relax, int depth, int start) : itsRelax(relax), itsDepth(depth), itsStart(start)
    { assert(depth>=2 && "PulayDMixer: a history shorter than 2 never extrapolates"); }

    //! \a working is this step's D_out and is NOT touched: the next Fock's D_in is mine (\c OwnedDIn).
    double Mix(const cd_t& working, const cd_t& old) override
    {
        const cd_t& in=InputFor(old);                                                // the D_in that drove THIS Fock
        const double dcd = working.GetChangeFrom(in)/working.GetTotalCharge();       // the legacy number (not the gate: D-MINDRHO)
        auto mixed=working.MakeCopy();
        if (++itsCount<=itsStart) { mixed->MixIn(in, 1.0-itsRelax); Install(std::move(mixed), working); return dcd; }   // priming
        if (itsStart>0 && itsCount==itsStart+1)                                      // R1.3: priming gives way to the history
            itsHandOff=HandOffRequest{"linear priming", "Pulay history on D (depth "+std::to_string(itsDepth)+")",
                                      "pulayStart "+std::to_string(itsStart)+" priming steps done"};
        itsIns.push_back(in.MakeCopy());                                             // D_in  of this step
        itsOuts.push_back(working.MakeCopy());                                       // D_out of this step
        while ((int)itsIns.size()>itsDepth) { itsIns.pop_front(); itsOuts.pop_front(); }
        const size_t n=itsIns.size();
        if (n<2) { mixed->MixIn(in, 1.0-itsRelax); Install(std::move(mixed), working); return dcd; }
        rsmat_t B=blazem::zero<double>(n);                                           // B_ij = <R_i,R_j>, bilinear in (in,out)
        for (size_t i=0;i<n;++i)
            for (size_t j=i;j<n;++j)
                B(i,j)= itsOuts[i]->Dot(*itsOuts[j]) - itsOuts[i]->Dot(*itsIns[j])
                      - itsIns [i]->Dot(*itsOuts[j]) + itsIns [i]->Dot(*itsIns[j]);
        const rvec_t c=qchem::Math::DIIS::Coefficients(qchem::Math::DIIS::Bordered(B));
        mixed->ReScale(0.0);                                                         // D_out is in the history
        for (size_t i=0;i<n;++i)
        {
            mixed->Axpy((1.0-itsRelax)*c[i], *itsIns [i]);
            mixed->Axpy(     itsRelax *c[i], *itsOuts[i]);
        }
        Install(std::move(mixed), working);
        return dcd;
    }
    const tChargeDensity<T>* FockDensity(const cd_t& working) const override { return &InputFor(working); }
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "PulD"; }
    bool HasGFilter() const override { return false; }   // mixes D itself
private:
    double itsRelax; int itsDepth, itsStart, itsCount=0;
    std::optional<HandOffRequest> itsHandOff;   //!< the priming -> history change, until the owner takes it (R1.3)
public:
    std::optional<HandOffRequest> TakeHandOff() override {auto h=std::move(itsHandOff); itsHandOff.reset(); return h;}
private:
    std::deque<std::unique_ptr<tMixableDensity<T>>> itsIns, itsOuts;   //!< the history (aligned index-wise)
};

} //namespace
