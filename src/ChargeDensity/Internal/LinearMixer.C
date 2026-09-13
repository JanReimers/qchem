// File: ChargeDensity/Internal/LinearMixer.C  Linear density-MATRIX mixing with the legacy adaptive α.
module;
#include <cmath>
export module qchem.ChargeDensity.Internal.LinearMixer;
export import qchem.ChargeDensity.DensityMixer;

export namespace qchem::ChargeDensity
{

//! Linear density-matrix mixing with the legacy adaptive-α heuristics.  α=1 (the default) = passthrough.
template <class T> class LinearMixer : public tDensityMixer<T>
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
    double GetRelax() const override { return itsRelax; }

    bool WantsReDamp(const MixSignals& s) const override
    {
        double dFD = s.FD - s.FDold;
        return s.FD>s.FDold && std::fabs(dFD)>1e-9;
    }
    double ReDampMix(cd_t& working, const cd_t& old) override
    {
        double dcd = working.GetChangeFrom(old);        // NOTE: un-normalised -- matches the legacy re-damp
        working.MixIn(old, 1.0-itsRelax/4.0);
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

} //namespace
