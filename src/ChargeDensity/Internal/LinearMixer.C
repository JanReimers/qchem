// File: ChargeDensity/Internal/LinearMixer.C  Linear density-MATRIX mixing with the legacy adaptive α.
module;
#include <cmath>
export module qchem.ChargeDensity.Internal.LinearMixer;
export import qchem.ChargeDensity.DensityMixer;

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

} //namespace
