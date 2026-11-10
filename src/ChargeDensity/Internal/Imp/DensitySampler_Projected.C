// File: ChargeDensity/Internal/Imp/DensitySampler_Projected.C  the PROJECTED strategy: XC sees P rho, H is <ij|P v>.
//
// V-XC-ONE-TERM step 2 (lane S, 2026-11-09): the plane-wave ball route, out of the pair sampler and into a strategy of its own.
// Until today a PW run built a PairDensitySampler and reached this route through two fallbacks: SampleOne's "the density cannot
// collocate, so take the ball round trip" and MatrixT's "this block has no raw adjoint, so fit v on the ball".  Both were
// decided per call from what the density or the block happened to answer.  Now the factory decides once, from what the fit
// basis states (its ball holds every orbital pair), and each strategy carries one route.  The arithmetic is the old ball
// route's, call for call, so every plane-wave number is unchanged.
module;
#include <mutex>
#include <cassert>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
module qchem.ChargeDensity.Internal.DensitySampler;
import qchem.Diagnostics;
import qchem.ChargeDensity;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.G_FieldEvaluator;
import qchem.Fitting.FunctionFitter;
import qchem.Reporting;
import qchem.Blaze;

namespace qchem::ChargeDensity
{

// Shared with the other two strategies (defined in the Singles / Pair units of this module).
bool HasExactSource(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd);
bool HasXCCorrection(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd);
void ReportNegativeRho(const DensitySamplerBase& q, const rvec_t& rho, const char* route);
bool BallHoldsPairs(const BasisSet::cFIT_SF_ABS& fit);

ProjectedDensitySampler::ProjectedDensitySampler(fbs_t fb)
    : itsFitBasis(std::move(fb))
    , itsExact(BallHoldsPairs(*itsFitBasis))
    , itsScalarFitter(Fitting::Factory(itsFitBasis))
    , itsRhoProjector(BasisSet::MakeGSpaceDensityProjector(itsFitBasis, "xc: rho-tilde IBZ star-average (SymmetrizeGMap)"))
{
    assert(itsFitBasis);
    if (!dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get()))
        throw std::logic_error("ProjectedDensitySampler: the fit basis has no raster (G_RasterTransform) to evaluate P rho on -- "
                               "only the plane-wave ball is realized on this strategy so far (V-XC-ONE-TERM step 2)");
    itsSource = itsExact ? XCRhoSource::ExactFit : XCRhoSource::Fitted;
}
ProjectedDensitySampler::~ProjectedDensitySampler() = default;

const BasisSet::G_RasterTransform& ProjectedDensitySampler::Raster() const
{
    return dynamic_cast<const BasisSet::G_RasterTransform&>(*itsFitBasis);   // checked at construction
}
double ProjectedDensitySampler::Integrate(const rvec_t& f) const {return Raster().Integral(f);}
size_t ProjectedDensitySampler::NumPoints() const {return Raster().RasterSize();}

// P rho of one density object.  Pin 28: a run that asked for the cusp deficit (XC from rho[D]) is served this route only where P
// is exact, since an exact P rho IS rho[D] and the deficit is identically zero (OpenWork §3, fixed 2026-10-26).
rvec_t ProjectedDensitySampler::SampleOne(const cChargeDensity* cd) const
{
    assert(cd);
    if (HasXCCorrection(cd) && !itsExact)
        throw std::logic_error("ProjectedDensitySampler: XCCuspDeficit asks for XC from rho[D] (doc/Pins.md pin 28), but this fit is "
            "inexact, so the projected route would feed XC the FITTED rho");
    return Raster().RhoOnGrid(cd->Project(*itsRhoProjector).GMap());
}

void ProjectedDensitySampler::Refresh(const cChargeDensity* cd) const
{
    assert(cd);
    assert(itsPolVersion==size_t(-1) && "ProjectedDensitySampler: this engine already served RhoPol -- the scalar and "
           "spin-resolved rho caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsRhoVersion) return;
    itsRhoVersion=cd->Version();
    itsRho=SampleOne(cd);
    if (qchem::Diagnostics::Enabled("xc_route"))
    {
        double rmin=1e300, rmax=-1e300;
        for (double r : itsRho) { rmin=std::min(rmin,r); rmax=std::max(rmax,r); }
        std::cout << "[xc route] PROJECTED (P rho on the ball, H = <ij|P v>; P " << (itsExact ? "exact" : "INEXACT") << ")"
                  << "  rho grid=[" << std::scientific << std::setprecision(2) << rmin << ", " << rmax << "]"
                  << " npts=" << itsRho.size() << std::defaultfloat << std::endl;
    }
    if (ReportGridCharge())
    {
        const double qGrid=Integrate(itsRho);
        const double qDM  =cd->GetTotalCharge();
        std::cout << "[grid charge] integral rho_grid=" << std::fixed << std::setprecision(6) << qGrid
                  << "  Tr(DS)=" << qDM
                  << "  lost=" << std::scientific << std::setprecision(3) << (qGrid-qDM)
                  << std::defaultfloat << std::endl;
    }
    ReportNegativeRho(*this, itsRho, "projected");
}

void ProjectedDensitySampler::RefreshPol(const cChargeDensity* cd) const
{
    assert(cd);
    assert(itsRhoVersion==size_t(-1) && "ProjectedDensitySampler: this engine already served the scalar Rho -- the two rho "
           "caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsPolVersion) return;
    itsPolVersion=cd->Version();
    const cChargeDensity* up=ChannelOf(cd, Spin::Up);
    const cChargeDensity* dn=ChannelOf(cd, Spin::Down);
    qchem::report::Timed timed("scf: XC projected rho sampling (per channel)");
    if (!up || !dn)
    {
        // The spin-agnostic seed: rho_up = rho_dn = rho/2 (the rule all three strategies apply).
        itsRhoUp=SampleOne(cd);
        itsRhoUp*=0.5;
        itsRhoDn=itsRhoUp;
        ReportNegativeRho(*this, itsRhoUp, "half-density seed");
        return;
    }
    // policy.xcFromDM (the wholesale retained-D repair) is wired only on the singles route; refused here as on the pair route.
    if ((HasExactSource(up) && !HasXCCorrection(up)) || (HasExactSource(dn) && !HasXCCorrection(dn)))
        throw std::logic_error("ProjectedDensitySampler: policy.xcFromDM is not wired on the projected route -- its rho comes from "
            "the density's projection onto the fit, not from a point projector, so the retained-D repair needs its own design.  Use "
            "the delta/singles quadrature for it.");
    itsRhoUp=SampleOne(up);
    itsRhoDn=SampleOne(dn);
    ReportNegativeRho(*this, itsRhoUp, "projected(up)");
    ReportNegativeRho(*this, itsRhoDn, "projected(dn)");
}

const rvec_t& ProjectedDensitySampler::Rho(const cChargeDensity* cd) const {Refresh(cd); return itsRho;}
const rvec_t& ProjectedDensitySampler::RhoPol(const cChargeDensity* cd, const Spin& s) const
{
    assert(s!=Spin::None && "ProjectedDensitySampler::RhoPol: ask for a channel, not the total");
    RefreshPol(cd);
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

chmat_t ProjectedDensitySampler::Matrix(const cobs_t* bs, const rvec_t& v) const
{
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(*bs);   // genuine "is it?" cross-cast (throws)
    std::lock_guard<std::mutex> lock(itsScalarFitterLock.m);   // ONE fitter refit per call: concurrent blocks take turns
    itsScalarFitter->DoFit(SampledField(v, NumPoints()));
    return itsScalarFitter->Overlap(orb);
}

rsmat_t ProjectedDensitySampler::Matrix(const robs_t*, const rvec_t&) const
{
    throw std::logic_error("ProjectedDensitySampler: a real TRIM block on the projected (ball-fit) XC route is not wired -- the ortho "
                           "scalar fitter has no FitContraction<double,dcmplx> face (V1.34)");
}

cvec_t ProjectedDensitySampler::Sample(const TransitionDensity<dcmplx>&, const Spin&) const
{
    throw std::logic_error("ProjectedDensitySampler::Sample: a linear response on the projected XC route is not realized -- a "
                           "transition density offers no projection onto the fit ball (V-XC-ONE-TERM)");
}
mat_t<dcmplx> ProjectedDensitySampler::Matrix(const cobs_t*, const cobs_t*, const cvec_t&, const rvec3_t&) const
{
    throw std::logic_error("ProjectedDensitySampler: the transition adjoint on the projected XC route is not realized (V-XC-ONE-TERM)");
}

} //namespace
