// File: ChargeDensity/Internal/Imp/DensitySampler_BallAtPoints.C  the plane-wave ball evaluated at a point set (I3).
//
// V-XC-ONE-TERM step 4 (lane S, 2026-11-14): a PW Vxc fit on a Becke XC mesh.  A decorator over the raster strategy (see the class note):
// its exact forward/adjoint pair, plus the exact pair raster -> ball -> points (forward FFT and truncation, then the direct synthesis)
// and back (the ball's projection from the points, then the inverse FFT).  First built ON the projected strategy's ball projector and
// fitter: exact on plane-wave orbitals (FD 7e-14) but 0.7% off dE/dD on a GPW block, whose ball projection and fitter contraction are
// not one pair there.  Composing with the raster strategy's own Rho/Matrix pair is what makes it exact on both.
module;
#include <cassert>
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

BallAtPointsDensitySampler::BallAtPointsDensitySampler(fbs_t fb, std::shared_ptr<const DensitySampler> raster, BasisSet::FitQuadrature points)
    : itsFitBasis(std::move(fb)), itsRaster(std::move(raster)), itsPoints(std::move(points))
{
    if (!itsFitBasis || !itsRaster) throw std::logic_error("BallAtPointsDensitySampler: needs the ball and its raster strategy");
    if (!itsPoints.GetMesh() || itsPoints.GetMesh()->size()==0) throw std::logic_error("BallAtPointsDensitySampler: a point set with no points");
    if (!dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get()) || !dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get()))
        throw std::logic_error("BallAtPointsDensitySampler: the fit basis must transform on a raster and evaluate / project its field at points");
}

double BallAtPointsDensitySampler::Integrate(const rvec_t& f) const
{
    const rvec_t& w=itsPoints.GetMesh()->Weights();
    assert(f.size()==w.size());
    double s=0.0;
    for (size_t g=0; g<w.size(); g++) s+=w[g]*f[g];
    return s;
}

// raster rho -> rho~ on the ball (the forward FFT, normalised /N, gathered on the ball's own {G}) -> the series at each point.
rvec_t BallAtPointsDensitySampler::AtPoints(const rvec_t& onRaster) const
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
rvec_t BallAtPointsDensitySampler::OnRaster(const rvec_t& atPoints) const
{
    qchem::report::Timed timed("scf: XC ball projection from the Becke points");
    const auto& raster=dynamic_cast<const BasisSet::G_RasterTransform&>(*itsFitBasis);
    const auto& field =dynamic_cast<const BasisSet::G_FieldEvaluator&>(*itsFitBasis);
    return raster.RhoOnGrid(field.ProjectField(itsPoints.GetMesh()->Points(), itsPoints.GetMesh()->Weights(), atPoints));
}

const rvec_t& BallAtPointsDensitySampler::Rho(const cChargeDensity* cd) const
{
    const rvec_t& r=itsRaster->Rho(cd);
    if (r.size()!=itsRhoFrom.size() || r!=itsRhoFrom) {itsRho=AtPoints(r); itsRhoFrom=r;}
    return itsRho;
}
const rvec_t& BallAtPointsDensitySampler::RhoPol(const cChargeDensity* cd, const Spin& s) const
{
    const rvec_t& r=itsRaster->RhoPol(cd, s);
    rvec_t& out =s==Spin::Up ? itsRhoUp : itsRhoDn;
    rvec_t& from=s==Spin::Up ? itsUpFrom : itsDnFrom;
    if (r.size()!=from.size() || r!=from) {out=AtPoints(r); from=r;}
    return out;
}

chmat_t BallAtPointsDensitySampler::Matrix(const cobs_t* bs, const rvec_t& v) const {return itsRaster->Matrix(bs, OnRaster(v));}
rsmat_t BallAtPointsDensitySampler::Matrix(const robs_t* bs, const rvec_t& v) const {return itsRaster->Matrix(bs, OnRaster(v));}

cvec_t BallAtPointsDensitySampler::Sample(const TransitionDensity<dcmplx>&, const Spin&) const
{
    throw std::logic_error("BallAtPointsDensitySampler::Sample: a linear response on the ball-at-Becke-points XC route is not realized (I3)");
}
mat_t<dcmplx> BallAtPointsDensitySampler::Matrix(const cobs_t*, const cobs_t*, const cvec_t&, const rvec3_t&) const
{
    throw std::logic_error("BallAtPointsDensitySampler: the transition adjoint on the ball-at-Becke-points XC route is not realized (I3)");
}

} //namespace
