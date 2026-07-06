// File: ChargeDensity/Imp/NumericCD.C  Superposition-of-atomic-densities DFT seed (impl).
module;
#include <cassert>
#include <vector>
#include <memory>
#include <cstddef>

module qchem.ChargeDensity.NumericCD;
import qchem.Blaze;                        // rvec3_t vector arithmetic (Gradient accumulation) + rvec_t
import qchem.BasisSet.Fit_IBS;             // FIT_SF_NonOrtho (the overlap metric-solve face)

namespace qchem::ChargeDensity
{

NumericCD::NumericCD(double totalCharge)
    : itsCharge(totalCharge), itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    assert(totalCharge>0);
}

void NumericCD::Insert(std::shared_ptr<const ScalarFunction<double>> d)
{
    assert(d);
    itsDensities.push_back(std::move(d));
    itsVersion=NextDensityVersion();   // mutated -> a new logical density
}

double NumericCD::operator()(const rvec3_t& r) const
{
    double rho=0;
    for (const auto& d : itsDensities) rho += (*d)(r);
    return itsScale*rho;
}

rvec3_t NumericCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (const auto& d : itsDensities) g += d->Gradient(r);
    return itsScale*g;
}

void NumericCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

// The seed's UNCONSTRAINED fit coefficients c0.  A numeric (SAD) seed has no density matrix, so its natural
// fit is the OVERLAP-metric fit of its own rho(r): c0 = S^-1 <f|rho> (sampling op(r)).  It is NOT a Coulomb
// fit -- so it overrides GetUnconstrainedFit directly (the ProjectedDensity_AO metric strategy) instead of
// faking a Coulomb RHS and letting the fitter's J^-1 cancel it.  Only the overlap face is needed: one honest
// "I want more" cross-cast to the S metric (item F removed the redundant J round-trip + Coulomb-face cast).
// The fitter then imposes the Dunlap charge constraint on top.
rvec_t NumericCD::GetUnconstrainedFit(const BasisSet::rFIT_CD_ABS* fbs) const
{
    const auto* sf = dynamic_cast<const BasisSet::FIT_SF_NonOrtho*>(fbs);   // the overlap metric-solve face
    assert(sf && "NumericCD::GetUnconstrainedFit: the seed overlap-fit needs the S metric (FIT_SF_NonOrtho face)");
    return sf->InvOverlap() * sf->Overlap(*this);   // c0 = S^-1 <f|rho>  (samples our rho(r))
}

} //namespace
