// File: ChargeDensity/Imp/NumericCD.C  Superposition-of-atomic-densities DFT seed (impl).
module;
#include <cassert>
#include <vector>
#include <memory>
#include <cstddef>

module qchem.ChargeDensity.NumericCD;
import qchem.Blaze;                        // rvec3_t vector arithmetic (Gradient accumulation) + rvec_t
import qchem.BasisSet.Fit_IBS;             // FIT_SF_NonOrtho / FIT_CD_NonOrtho (the two metric-solve faces)

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

// The seed's own density-fit projection <rho|c>.  A numeric (SAD) seed has no density matrix, so we overlap-
// fit our own rho(r) (e = S^-1 <f|rho>, sampling op(r)) and Coulomb-project (J e).  The downstream
// ConstrainedFF then solves c0 = J^-1 (J e) = e + the Dunlap charge constraint -- bit-identical to the old
// FittedCD seed path.  (Relocated from FittedCDImp::ScalarSeedProjection_AO; still carries the CD->SF cross-
// cast to reach the overlap metric -- that is item F in doc/FittingCleanupPlan.md.)
rvec_t NumericCD::GetRepulsion3C(const BasisSet::rFIT_CD_ABS* fbs) const
{
    const auto* sf = dynamic_cast<const BasisSet::FIT_SF_NonOrtho*>(fbs);   // the overlap metric-solve face
    const auto* no = dynamic_cast<const BasisSet::FIT_CD_NonOrtho*>(fbs);   // the Coulomb metric face
    assert(sf && "NumericCD::GetRepulsion3C: the CD fit basis must also expose its overlap-metric (FIT_SF_NonOrtho) face");
    assert(no && "NumericCD::GetRepulsion3C: the seed overlap-fit needs the Coulomb metric (FIT_CD_NonOrtho)");
    rvec_t e = sf->InvOverlap() * sf->Overlap(*this);   // overlap-metric fit coeffs (samples our rho(r))
    return no->Repulsion() * e;                         // Coulomb-project -> <rho_fit|f_c>
}

} //namespace
