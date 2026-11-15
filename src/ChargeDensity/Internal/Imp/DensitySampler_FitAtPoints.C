// File: ChargeDensity/Internal/Imp/DensitySampler_FitAtPoints.C  XC of a fitted density at a point set: the decorator and its fit spaces.
//
// V-XC-ONE-TERM step 4 (lane S, 2026-11-14): a PW Vxc fit on a Becke XC mesh.  A decorator over the raster strategy (see the class note):
// its exact forward/adjoint pair, plus the exact pair raster -> ball -> points (forward FFT and truncation, then the direct synthesis)
// and back (the ball's projection from the points, then the inverse FFT).  First built ON the projected strategy's ball projector and
// fitter: exact on plane-wave orbitals (FD 7e-14) but 0.7% off dE/dD on a GPW block, whose ball projection and fitter contraction are
// not one pair there.  Composing with the raster strategy's own Rho/Matrix pair is what makes it exact on both.
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
module qchem.ChargeDensity.Internal.DensitySampler;
import qchem.ChargeDensity;
import qchem.BasisSet.G_FieldEvaluator;
import qchem.Reporting;
import qchem.Parallel;    // InnerThreads: the synthesis is per-point parallel
import qchem.Blaze;

namespace qchem::ChargeDensity
{

// ---- the plane-wave ball as a fit space ---------------------------------------------------------------------------------------------
BallFitSpace::BallFitSpace(std::shared_ptr<const BasisSet::cFIT_SF_ABS> fb, BasisSet::FitQuadrature points)
    : itsFitBasis(std::move(fb)), itsPoints(std::move(points))
{
    if (!itsFitBasis) throw std::logic_error("BallFitSpace: needs the ball");
    if (!itsPoints.GetMesh() || itsPoints.GetMesh()->size()==0) throw std::logic_error("BallFitSpace: a point set with no points");
    if (!dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get()) || !dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get()))
        throw std::logic_error("BallFitSpace: the fit basis must transform on a raster and evaluate / project its field at points");
}

// ---- the decorator --------------------------------------------------------------------------------------------------------------------
FitAtPointsDensitySampler::FitAtPointsDensitySampler(std::shared_ptr<const FitSpace> fit, std::shared_ptr<const DensitySampler> raster, BasisSet::FitQuadrature points)
    : itsFit(std::move(fit)), itsRaster(std::move(raster)), itsPoints(std::move(points))
{
    if (!itsFit || !itsRaster) throw std::logic_error("FitAtPointsDensitySampler: needs the fit space and its raster strategy");
    if (!itsPoints.GetMesh() || itsPoints.GetMesh()->size()==0) throw std::logic_error("FitAtPointsDensitySampler: a point set with no points");
}

double FitAtPointsDensitySampler::Integrate(const rvec_t& f) const
{
    const rvec_t& w=itsPoints.GetMesh()->Weights();
    assert(f.size()==w.size());
    double s=0.0;
    for (size_t g=0; g<w.size(); g++) s+=w[g]*f[g];
    return s;
}

// raster rho -> rho~ on the ball (the forward FFT, normalised /N, gathered on the ball's own {G}) -> the series at each point.
rvec_t BallFitSpace::AtPoints(const rvec_t& onRaster) const
{
    qchem::report::Timed timed("scf: XC ball synthesis at the Becke points");
    const auto& raster=dynamic_cast<const BasisSet::G_RasterTransform&>(*itsFitBasis);   // checked at construction
    const auto& field =dynamic_cast<const BasisSet::G_FieldEvaluator&>(*itsFitBasis);
    const ΔG_Map rhoTilde=raster.FieldCoeffs(raster.ForwardFFT(onRaster));
    const rvec3vec_t& R=itsPoints.GetMesh()->Points();
    rvec_t rho(R.size());
#ifdef QCHEM_OPENMP
    #pragma omp parallel for schedule(static) num_threads(qchem::InnerThreads())
#endif
    for (size_t g=0; g<R.size(); g++) rho[g]=field.EvalField(rhoTilde, R[g]);
    return rho;
}

// v at my points -> v~(G) = (1/Omega) sum_g w_g v_g e^{-iG.r_g} on the ball -> the band-limited field on the raster.  The raster strategy's
// Matrix(v_raster) is sum_k (Omega/N) v_k d rho_k/dD, and dE/d rho_k = (Omega/N) v_fit(r_k) for E = sum_g w_g e(rho_P(r_g)): no other factor.
rvec_t BallFitSpace::OnRaster(const rvec_t& atPoints) const
{
    qchem::report::Timed timed("scf: XC ball projection from the Becke points");
    const auto& raster=dynamic_cast<const BasisSet::G_RasterTransform&>(*itsFitBasis);
    const auto& field =dynamic_cast<const BasisSet::G_FieldEvaluator&>(*itsFitBasis);
    return raster.RhoOnGrid(field.ProjectField(itsPoints.GetMesh()->Points(), itsPoints.GetMesh()->Weights(), atPoints));
}

void FitAtPointsDensitySampler::ReportHealth(const rvec_t& rhoP, double N) const
{
    double minRho=rhoP.size() ? rhoP[0] : 0.0;
    for (double x : rhoP) minRho=std::min(minRho, x);
    const double dq=Integrate(rhoP)-N;
    const SamplingTarget* t=itsPoints.GetMesh()->FourierTarget();
    const double target= t ? t->accuracy.target : 0.0;   // 0 = the run stated no target: report, and warn on a negative density only
    {
        qchem::report::Scope here("xcFitAtPoints");
        qchem::report::Set("fitSpace", itsFit->Name());
        qchem::report::Set("chargeError", dq);
        qchem::report::Set("minRho", minRho);
        qchem::report::Set("chargeTarget", target);
    }
    const bool overCharge= target>0.0 && std::fabs(dq)>target, negative= minRho<0.0;
    if (!itsAnnounced || (overCharge && !itsWarnedCharge) || (negative && !itsWarnedNegative))
    {
        std::cout << "[xc fit at points] " << itsFit->Name() << ": charge error " << dq << " e (target " << target << ")  min rho_P " << minRho;
        if (overCharge) std::cout << "  WARNING: the fit misses the charge target (no eps sizing of the fit yet)";
        if (negative)   std::cout << "  WARNING: rho_P < 0 at some points (a sharp cut rings, pin 1; a fit undershoots): XC drops them";
        std::cout << std::endl;
        itsAnnounced=true; itsWarnedCharge|=overCharge; itsWarnedNegative|=negative;
    }
}

const rvec_t& FitAtPointsDensitySampler::Rho(const cChargeDensity* cd) const
{
    const rvec_t& r=itsRaster->Rho(cd);
    if (r.size()!=itsRhoFrom.size() || r!=itsRhoFrom) {itsRho=itsFit->AtPoints(r); itsRhoFrom=r; ReportHealth(itsRho, cd->GetTotalCharge());}
    return itsRho;
}
const rvec_t& FitAtPointsDensitySampler::RhoPol(const cChargeDensity* cd, const Spin& s) const
{
    const rvec_t& r=itsRaster->RhoPol(cd, s);
    rvec_t& out =s==Spin::Up ? itsRhoUp : itsRhoDn;
    rvec_t& from=s==Spin::Up ? itsUpFrom : itsDnFrom;
    if (r.size()!=from.size() || r!=from)
    {
        out=itsFit->AtPoints(r); from=r;
        // the health of the total: once both channels are synthesized for the same density
        if (s==Spin::Down && itsRhoUp.size()==out.size()) ReportHealth(rvec_t(itsRhoUp+out), cd->GetTotalCharge());
    }
    return out;
}

chmat_t FitAtPointsDensitySampler::Matrix(const cobs_t* bs, const rvec_t& v) const {return itsRaster->Matrix(bs, itsFit->OnRaster(v));}
rsmat_t FitAtPointsDensitySampler::Matrix(const robs_t* bs, const rvec_t& v) const {return itsRaster->Matrix(bs, itsFit->OnRaster(v));}

cvec_t FitAtPointsDensitySampler::Sample(const TransitionDensity<dcmplx>&, const Spin&) const
{
    throw std::logic_error("FitAtPointsDensitySampler::Sample: a linear response on the fit-at-points XC route is not realized");
}
mat_t<dcmplx> FitAtPointsDensitySampler::Matrix(const cobs_t*, const cobs_t*, const cvec_t&, const rvec3_t&) const
{
    throw std::logic_error("FitAtPointsDensitySampler: the transition adjoint on the fit-at-points XC route is not realized");
}

} //namespace
