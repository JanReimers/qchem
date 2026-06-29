// File: ChargeDensity/Imp/CompositeFittedCD.C  Superposition-of-atomic-densities DFT seed (impl).
module;
#include <cassert>
#include <vector>
#include <memory>
#include <cstddef>

module qchem.ChargeDensity.CompositeFittedCD;
import qchem.BasisSet.Fit_IBS;             // FIT_SF_ABS (overlap face of the term's fit basis)
import qchem.Blaze;                        // rvec_t, rsmat_t, matrix-vector products

namespace qchem::ChargeDensity
{

CompositeFittedCD::CompositeFittedCD(double totalCharge)
    : itsCharge(totalCharge), itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    assert(totalCharge>0);
}

void CompositeFittedCD::Insert(std::shared_ptr<const ScalarFunction<double>> d)
{
    assert(d);
    itsDensities.push_back(std::move(d));
    itsVersion=NextDensityVersion();   // mutated -> a new logical density
}

double CompositeFittedCD::operator()(const rvec3_t& r) const
{
    double rho=0;
    for (const auto& d : itsDensities) rho += (*d)(r);
    return itsScale*rho;
}

rvec3_t CompositeFittedCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (const auto& d : itsDensities) g += d->Gradient(r);
    return itsScale*g;
}

void CompositeFittedCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

// <rho|c> = the Coulomb projection of this real-space density onto fit basis \a fbs.  Derived on demand
// from op(r): overlap-fit this density onto fbs (e = S^-1 <f|rho>, the FIT_SF_ABS face), then Coulomb-
// project (Repulsion * e, the FIT_CD_ABS face).  fbs is concretely a Fit_IBS, so the FIT_CD_ABS -> FIT_SF_ABS
// cross-cast (abstract-to-abstract) succeeds.  FittedVee then recovers e and builds the Hartree matrix.
rvec_t CompositeFittedCD::GetRepulsion3C(const BasisSet::FIT_CD_ABS* fbs) const
{
    const auto* sf = dynamic_cast<const BasisSet::FIT_SF_ABS*>(fbs);
    assert(sf && "CompositeFittedCD: the CD fit basis must also expose its overlap (FIT_SF_ABS) face");
    rvec_t e = sf->InvOverlap() * sf->Overlap(*this);   // overlap-metric fit coeffs (uses op(r), incl. itsScale)
    return fbs->Repulsion() * e;                         // <rho_fit|f_c>, Coulomb metric
}

} //namespace
