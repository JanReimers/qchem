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
#include <vector>
module qchem.ChargeDensity.Internal.DensitySampler;
import qchem.ChargeDensity;
import qchem.BasisSet.G_FieldEvaluator;
import qchem.BasisSet.Orbital_1E_IBS;   // the periodic A1 Gamma block's overlap and point values
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
        itsFit->ReportHealthFacts();
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

// ---- the periodic A1 Gaussian set as a fit space (V-XC-ONE-TERM step 6) ----------------------------------------------------------------
namespace qchem::ChargeDensity
{

GaussianFitSpace::GaussianFitSpace(std::shared_ptr<const BasisSet::Orbital_1E_IBS<double>> a1, BasisSet::FitQuadrature points)
    : itsA1(std::move(a1)), itsPoints(std::move(points))
{
    if (!itsA1) throw std::logic_error("GaussianFitSpace: needs the periodic A1 block");
    if (!itsPoints.GetMesh() || itsPoints.GetMesh()->size()==0) throw std::logic_error("GaussianFitSpace: a point set with no points");
    // S^-1, robust to a near-dependent set: S = V diag(w) V^T, the modes below kMetricRelTol x w_max dropped and REPORTED.
    const auto& S=itsA1->Overlap();
    const size_t n=S.rows();
    rsmat_t Sd(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) Sd(i,j)=S(i,j);
    rvec_t w; rmat_t V;
    blazem::eigen(Sd, w, V);
    double wmax=0.0;
    for (size_t k=0;k<n;k++) wmax=std::max(wmax, w[k]);
    size_t dropped=0;
    for (size_t k=0;k<n;k++) if (w[k]<kMetricRelTol*wmax) dropped++;
    if (dropped)
    {
        qchem::report::Set("a1MetricModesDropped", (long)dropped);
        qchem::report::Set("a1MetricFloor", kMetricRelTol*wmax);
        std::cout<<"[periodic A1 fit] the overlap metric of "<<n<<" functions is near-singular: "<<dropped<<" mode(s) below the floor "
                 <<kMetricRelTol*wmax<<" (= "<<kMetricRelTol<<" x the largest eigenvalue) are dropped"<<std::endl;
    }
    itsSinv=rsmat_t(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            double s=0.0;
            for (size_t k=0;k<n;k++) if (w[k]>=kMetricRelTol*wmax) s+=V(i,k)*V(j,k)/w[k];
            itsSinv(i,j)=s;
        }
    // THE METRIC GAP: the XC quadrature's metric Phi^T w Phi against the analytic S, streamed in chunks.
    {
        const rvec3vec_t& R=itsPoints.GetMesh()->Points();
        const rvec_t&     W=itsPoints.GetMesh()->Weights();
        rmat_t Sq(n, n, 0.0);
        for (size_t lo=0; lo<R.size(); lo+=2048)
        {
            const size_t hi=std::min(R.size(), lo+2048);
            rvec3vec_t sub(hi-lo);
            for (size_t g=lo; g<hi; g++) sub[g-lo]=R[g];
            const rmat_t P=(*itsA1)(sub);
            for (size_t g=lo; g<hi; g++) for (size_t a=0; a<n; a++) { const double wa=W[g]*P(g-lo,a); for (size_t b=0; b<n; b++) Sq(a,b)+=wa*P(g-lo,b); }
        }
        double num=0.0, den=0.0;
        for (size_t a=0; a<n; a++) for (size_t b=0; b<n; b++) { const double d=S(a,b)-Sq(a,b); num+=d*d; den+=S(a,b)*S(a,b); }
        itsMetricGap=std::sqrt(num/std::max(den, 1e-300));
        qchem::report::Set("a1MetricGap", itsMetricGap);
        std::cout<<"[periodic A1 fit] metric gap ||S - Phi^T w Phi|| / ||S|| = "<<itsMetricGap<<" on "<<R.size()<<" XC points"
                 <<(itsMetricGap>1e-3 ? "  WARNING: the XC mesh does not resolve the set's functions (R is not a projector there)" : "")<<std::endl;
    }
    // The footprint facts, one level in: what streaming saves, against what a stored table would cost (the policy layer judges them).
    qchem::report::Set("a1Functions", (long)n);
    qchem::report::Set("a1XCPoints", (long)itsPoints.GetMesh()->size());
}

// Phi(R) evaluated in chunks of points (the block's image-summed bulk evaluation): never the whole table at once.
static constexpr size_t kA1Chunk=2048;
rvec_t GaussianFitSpace::PhiT(const rvec3vec_t& R, const rvec_t& f) const
{
    assert(f.size()==R.size());
    const size_t n=NumFunctions(), nblk=(R.size()+kA1Chunk-1)/kA1Chunk;
    std::vector<rvec_t> part(nblk, rvec_t(n, 0.0));   // one partial per chunk: a deterministic sum at any thread count
#ifdef QCHEM_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(qchem::InnerThreads())
#endif
    for (size_t b=0; b<nblk; b++)
    {
        const size_t lo=b*kA1Chunk, hi=std::min(R.size(), lo+kA1Chunk);
        rvec3vec_t sub(hi-lo);
        for (size_t g=lo; g<hi; g++) sub[g-lo]=R[g];
        const rmat_t P=(*itsA1)(sub);
        for (size_t g=lo; g<hi; g++) for (size_t a=0; a<n; a++) part[b][a]+=P(g-lo,a)*f[g];
    }
    rvec_t out(n, 0.0);
    for (const auto& p : part) out+=p;
    return out;
}
rvec_t GaussianFitSpace::Phi(const rvec3vec_t& R, const rvec_t& c) const
{
    const size_t n=NumFunctions(), nblk=(R.size()+kA1Chunk-1)/kA1Chunk;
    rvec_t out(R.size(), 0.0);
#ifdef QCHEM_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(qchem::InnerThreads())
#endif
    for (size_t b=0; b<nblk; b++)
    {
        const size_t lo=b*kA1Chunk, hi=std::min(R.size(), lo+kA1Chunk);
        rvec3vec_t sub(hi-lo);
        for (size_t g=lo; g<hi; g++) sub[g-lo]=R[g];
        const rmat_t P=(*itsA1)(sub);
        for (size_t g=lo; g<hi; g++) { double s=0.0; for (size_t a=0; a<n; a++) s+=P(g-lo,a)*c[a]; out[g]=s; }
    }
    return out;
}

rvec_t GaussianFitSpace::Resample(const rvec_t& atPoints) const
{
    qchem::report::Timed timed("scf: XC periodic A1 fit of v");
    const rvec_t& w=itsPoints.GetMesh()->Weights();
    rvec_t wv(atPoints.size());
    for (size_t g=0; g<wv.size(); g++) wv[g]=w[g]*atPoints[g];
    const rvec_t c=itsSinv*PhiT(itsPoints.GetMesh()->Points(), wv);
    return Phi(itsPoints.GetMesh()->Points(), c);
}

void GaussianFitSpace::ReportHealthFacts() const {qchem::report::Set("metricGap", itsMetricGap);}

// ---- the fitted-potential cell -------------------------------------------------------------------------------------------------------
FittedPotentialDensitySampler::FittedPotentialDensitySampler(std::shared_ptr<const GaussianFitSpace> fit, std::shared_ptr<const DensitySampler> inner)
    : itsFit(std::move(fit)), itsInner(std::move(inner))
{
    if (!itsFit || !itsInner) throw std::logic_error("FittedPotentialDensitySampler: needs the fit space and its inner strategy");
}
cvec_t FittedPotentialDensitySampler::Sample(const TransitionDensity<dcmplx>&, const Spin&) const
{
    throw std::logic_error("FittedPotentialDensitySampler::Sample: a linear response on the fitted-potential cell is not realized");
}
mat_t<dcmplx> FittedPotentialDensitySampler::Matrix(const cobs_t*, const cobs_t*, const cvec_t&, const rvec3_t&) const
{
    throw std::logic_error("FittedPotentialDensitySampler: the transition adjoint on the fitted-potential cell is not realized");
}

} //namespace
