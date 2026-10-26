// File: ChargeDensity/Internal/LinearMixer.C  The density-MATRIX mixers: linear (legacy adaptive α) and D-space Pulay.
module;
#include <cmath>
#include <deque>
#include <memory>
#include <cassert>
export module qchem.ChargeDensity.Internal.LinearMixer;
export import qchem.ChargeDensity.DensityMixer;
import qchem.Math.DIIS;   // Bordered / Coefficients -- the shared extrapolator (the ρ̃ Pulay and the Fock DIIS use it too)
import qchem.Blaze;

export namespace qchem::ChargeDensity
{

//! Linear density-matrix mixing with the legacy adaptive-α heuristics.  α=1 (the default) = passthrough.
template <class T> class LinearMixer : public tDensityMixer<T>, public virtual tAdaptiveMixer<T>
{
public:
    typedef typename tDensityMixer<T>::cd_t cd_t;
    explicit LinearMixer(double relax0) : itsRelax(relax0) {}

    double Mix(cd_t& working, const cd_t& old) override
    {
        double dcd = working.GetChangeFrom(old)/working.GetTotalCharge();      // relative MaxAbs change
        if (dcd<1e-5) itsRelMax=0.5;
        working.MixIn(old, 1.0-itsRelax);                                      // (1−relax)ρ_in + relax ρ_out
        return dcd;
    }
    const tChargeDensity<T>* FockDensity(const cd_t& working) const override { return &working; }   // D-mixing: the mixed D IS the Fock density
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "Lin"; }
    bool HasGFilter() const override { return false; }   // mixes D itself

    //! The legacy [F,D]-keyed policy, in one place: if this step made [F,D] WORSE, retract it to a quarter step
    //! and shrink α; if it improved, grow α (clamped).
    //!
    //! THE RETRACTION NEEDS NO REBUILT ρ_out.  \a working already holds ρ_w = (1−α)ρ_in + αρ_out, and the mix
    //! is linear, so the quarter step (1−α')ρ_in + α'ρ_out with α' = α/4 is
    //!     ρ_in + (α'/α)(ρ_w − ρ_in)  =  MixIn(old, 1 − α'/α)  =  MixIn(old, 0.75),
    //! and the legacy gate ‖ρ_out − ρ_in‖ (un-normalised, as it always was) is ‖ρ_w − ρ_in‖/α.  The old
    //! choreography had the iterator rebuild ρ_out from the wave function for this -- a second full density
    //! build -- and the numbers differ from that only at rounding level.
    bool Adapt(const MixSignals& s, cd_t& working, const cd_t& old, double& dRho) override
    {
        bool remixed=false;
        const double dFD = s.FD - s.FDold;
        if (s.FD>s.FDold && std::fabs(dFD)>1e-9 && itsRelax>0.0)
        {
            dRho = working.GetChangeFrom(old)/itsRelax;
            working.MixIn(old, 0.75);
            itsRelax*=0.8;
            remixed=true;
        }
        if (s.FD<s.FDold) itsRelax*=1.5;
        if (itsRelax>itsRelMax) itsRelax=itsRelMax;
        return remixed;
    }
private:
    double itsRelax;
    double itsRelMax=1.0;
};

//! \brief PULAY (density-DIIS) ON THE DENSITY MATRIX -- "Pulay without Kerker" (OpenWork §3, 2026-10-26).
//!
//! The history is D itself: copies of each step's D_in and D_out.  The step solves the DIIS bordered system over the
//! residuals R_i = D_out,i - D_in,i (inner product \c tMixableDensity::Dot) for c (Σc = 1), and the next input is
//! \f$D_{in}^{next}=(1-\alpha)\sum_i c_iD_{in,i}+\alpha\sum_i c_iD_{out,i}\f$, formed IN PLACE on the working density with
//! \c ReScale + \c Axpy.  So the Fock -- and V_xc -- see that mixed density MATRIX through the DM route: no fit, no ρ̃, no
//! raster shadow, no cusp deficit (pin 28 holds by construction; the loop's RequirePin28 would throw otherwise).  ρ[D] is
//! linear in D, so this IS Pulay on ρ[D] -- exactly, with no G-ball truncation.  Any T, any system: it asks nothing of the
//! geometry.  On a polarized run the composite holds both spins, so ONE history and one B summed over them come free.
//! Priming: \a start plain linear steps at \a relax before the history engages (history-based mixing is unstable far out).
//! Memory: 2 copies of D per history slot (the residual inner products are formed bilinearly from in/out).
template <class T> class PulayDMixer : public tDensityMixer<T>
{
public:
    typedef typename tDensityMixer<T>::cd_t cd_t;
    PulayDMixer(double relax, int depth, int start) : itsRelax(relax), itsDepth(depth), itsStart(start)
    { assert(depth>=2 && "PulayDMixer: a history shorter than 2 never extrapolates"); }

    double Mix(cd_t& working, const cd_t& old) override
    {
        const double dcd = working.GetChangeFrom(old)/working.GetTotalCharge();   // the legacy number (not the gate: D-MINDRHO)
        if (++itsCount<=itsStart) { working.MixIn(old, 1.0-itsRelax); return dcd; }  // priming: the linear step
        itsIns.push_back(old.MakeCopy());                                            // D_in  of this step
        itsOuts.push_back(working.MakeCopy());                                       // D_out of this step
        while ((int)itsIns.size()>itsDepth) { itsIns.pop_front(); itsOuts.pop_front(); }
        const size_t n=itsIns.size();
        if (n<2) { working.MixIn(old, 1.0-itsRelax); return dcd; }
        rsmat_t B=blazem::zero<double>(n);                                           // B_ij = <R_i,R_j>, bilinear in (in,out)
        for (size_t i=0;i<n;++i)
            for (size_t j=i;j<n;++j)
                B(i,j)= itsOuts[i]->Dot(*itsOuts[j]) - itsOuts[i]->Dot(*itsIns[j])
                      - itsIns [i]->Dot(*itsOuts[j]) + itsIns [i]->Dot(*itsIns[j]);
        const rvec_t c=qchem::Math::DIIS::Coefficients(qchem::Math::DIIS::Bordered(B));
        working.ReScale(0.0);                                                        // D_out is in the history now
        for (size_t i=0;i<n;++i)
        {
            working.Axpy((1.0-itsRelax)*c[i], *itsIns [i]);
            working.Axpy(     itsRelax *c[i], *itsOuts[i]);
        }
        return dcd;
    }
    const tChargeDensity<T>* FockDensity(const cd_t& working) const override { return &working; }   // the mixed D IS the Fock density
    double GetRelax() const override { return itsRelax; }
    const char* Tag() const override { return "PulD"; }
    bool HasGFilter() const override { return false; }   // mixes D itself
private:
    double itsRelax; int itsDepth, itsStart, itsCount=0;
    std::deque<std::unique_ptr<tMixableDensity<T>>> itsIns, itsOuts;   //!< the history (aligned index-wise)
};

} //namespace
