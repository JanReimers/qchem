// File: ChargeDensity/Imp/CompositeFittedCD.C  Superposition-of-atomic-densities DFT seed (impl).
module;
#include <cassert>
#include <vector>
#include <memory>
#include <cstddef>

module qchem.ChargeDensity.CompositeFittedCD;
import qchem.ChargeDensity.Types;          // ohfbs_t
import qchem.BasisSet.Fit_IBS;             // FIT_SF_ABS (overlap face of the term's fit basis)
import qchem.Blaze;                        // rvec_t, rsmat_t, matrix-vector products

namespace qchem::ChargeDensity
{

static size_t NextSeedVersion() {static size_t n=0; return ++n;}   // transient freshness serial

CompositeFittedCD::CompositeFittedCD(double totalCharge)
    : itsCharge(totalCharge), itsVersion(NextSeedVersion())
{
    assert(totalCharge>0);
}

void CompositeFittedCD::Insert(std::shared_ptr<const ScalarFunction<double>> d)
{
    assert(d);
    itsDensities.push_back(std::move(d));
    itsVersion=NextSeedVersion();   // mutated -> a new logical density
}

double CompositeFittedCD::operator()(const rvec3_t& r) const
{
    double rho=0;
    for (const auto& d : itsDensities) rho += (*d)(r);
    return rho;
}

rvec3_t CompositeFittedCD::Gradient(const rvec3_t& r) const
{
    rvec3_t g(0,0,0);
    for (const auto& d : itsDensities) g += d->Gradient(r);
    return g;
}

// <rho|c> = the Coulomb projection of this real-space density onto fit basis \a fbs.  Derived on demand
// from op(r): overlap-fit this density onto fbs (e = S^-1 <f|rho>, the FIT_SF_ABS face), then Coulomb-
// project (Repulsion * e, the FIT_CD_ABS face).  fbs is concretely a Fit_IBS, so the FIT_CD_ABS -> FIT_SF_ABS
// cross-cast (abstract-to-abstract) succeeds.  FittedVee then recovers e and builds the Hartree matrix.
rvec_t CompositeFittedCD::GetRepulsion3C(const BasisSet::FIT_CD_ABS* fbs) const
{
    const auto* sf = dynamic_cast<const BasisSet::FIT_SF_ABS*>(fbs);
    assert(sf && "CompositeFittedCD: the CD fit basis must also expose its overlap (FIT_SF_ABS) face");
    rvec_t e = sf->InvOverlap() * sf->Overlap(*this);   // overlap-metric fit coeffs (uses op(r))
    return fbs->Repulsion() * e;                         // <rho_fit|f_c>, Coulomb metric
}

//---- density-MATRIX capabilities: a fitted seed has none (deferred tChargeDensity/tDM_CD ISP split) ----
// These are unreachable on a seed: it is consumed once, at iteration 0, by the DFT Fock build (op(r) for
// Vxc, GetRepulsion3C for Vee) -- never mixed, change-measured, energy-contracted, or used for HF J/K.

double CompositeFittedCD::DM_Contract(const tStatic_CC<double>*) const
{
    assert(false && "CompositeFittedCD::DM_Contract -- DFT seed has no density matrix");
    return 0.0;
}
double CompositeFittedCD::DM_Contract(const tDynamic_CC<double>*,const tDM_CD<double>*) const
{
    assert(false && "CompositeFittedCD::DM_Contract -- DFT seed has no density matrix");
    return 0.0;
}
void CompositeFittedCD::ReScale(double)
{
    assert(false && "CompositeFittedCD::ReScale -- DFT seed is not mixed/rescaled in the SCF loop");
}
void CompositeFittedCD::MixIn(const tDM_CD<double>&,double)
{
    assert(false && "CompositeFittedCD::MixIn -- DFT seed is not mixed in the SCF loop");
}
double CompositeFittedCD::GetChangeFrom(const tDM_CD<double>&) const
{
    assert(false && "CompositeFittedCD::GetChangeFrom -- DFT seed is not convergence-checked");
    return 0.0;
}
void CompositeFittedCD::AccumulateDirect(hmat_t<double>&, const ohfbs_t*) const
{
    assert(false && "CompositeFittedCD::AccumulateDirect -- HF J needs a density matrix; SAD seed is DFT-only");
}
void CompositeFittedCD::AccumulateExchange(hmat_t<double>&, const ohfbs_t*) const
{
    assert(false && "CompositeFittedCD::AccumulateExchange -- HF K needs a density matrix; SAD seed is DFT-only");
}

} //namespace
