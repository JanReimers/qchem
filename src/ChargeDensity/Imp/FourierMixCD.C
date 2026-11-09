// File: ChargeDensity/Imp/FourierMixCD.C  Kerker-mixable G-space density (see FourierMixCD.C).
module;
#include <algorithm>   // std::min/max (the batched op()'s dm ranges)
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>     // getenv (GPW_KERKER_SPECTRUM)
#include <iomanip>
#include <iostream>    // the residual-spectrum census goes to cout, like the other GPW_* instruments
#include <map>
#include <memory>      // shared_ptr (the deposited cusp-deficit correction field)
#include <stdexcept>
#include <vector>
module qchem.ChargeDensity.FourierMixCD;
import qchem.Types;        // rvec3_t, dcmplx
import qchem.PointSampling; // FourierSampler, MakeNUFFTSampler (the batched sampler)
import qchem.Parallel;     // WorkerThreads (QCHEM_OPENMP_THREADS -- the batched inverse FT over mesh points)

namespace qchem::ChargeDensity
{

FourierMixCD::FourierMixCD(ΔG_Map rhoTilde, ReciprocalLattice recip, double charge)
    : itsRho(std::move(rhoTilde)), itsRecip(std::move(recip))
    , itsCharge(charge), itsVersion(NextDensityVersion())
{}

FourierMixCD::FourierMixCD(ΔG_Map rhoTilde, ReciprocalLattice recip, double charge, Extras extras)
    : itsRho(std::move(rhoTilde)), itsRhoRaw(std::move(extras.raster))
    , itsXCCorrection(std::move(extras.xcCorrection)), itsEffectiveAlpha(extras.alphaEff)
    , itsRecip(std::move(recip)), itsCharge(charge), itsVersion(NextDensityVersion())
{}

ΔG_Map FourierMixCD::FourierComponents(const BasisSet::cFIT_SF_ABS&) const
{
    if (itsScale==1.0) return itsRho;
    ΔG_Map r=itsRho; for (auto& [dm,v] : r) v = itsScale*dcmplx(v); return r;   // apply ReScale on the way out
}

double FourierMixCD::operator()(const rvec3_t& r) const
{
    // rho(r) = Sum_G rho-tilde(G) e^{iG.r} (inverse FT).  Real by construction for a real density
    // (conjugate-symmetric rho-tilde).  Single-point form; the mesh-XC quadrature batches via the
    // rvec3vec_t overload below.
    dcmplx s(0.0);
    for (const auto& [dm, rt] : itsRho)
    {
        rvec3_t G = itsRecip.GetCell().ToCartesian(rvec3_t(dm));
        s += dcmplx(rt) * std::exp(dcmplx(0.0, G*r));
    }
    return itsScale*std::real(s);
}

// THE rho-tilde SAMPLER (doc/OpenWork.md §4, step 1, user 2026-11-03/04).  rho(r) = Sum_dm c(dm) e^{iG.r} with G = B dm is
// evaluated by a sampler the POINT SET makes (qchem.PointSampling): a type-2 NUFFT on any set, the exact inverse FFT on
// the uniform cell raster.  No lattice sum is cut (pin 1).
//
// WHY: the direct per-point G-sum this replaces was O(N_G * N_r) and grew with the density ball though the mesh does not
// -- the matrix-free and the N4 cusp-deficit samplings, both through here, were 1415 + 689 s of a 4100 s MnO geometry
// SCF.  The pointwise operator() above is the direct series, kept as the reference.
static FourierBox BoxOf(const ΔG_Map& rho, const ReciprocalLattice& recip)
{
    ivec3_t lo=rho.begin()->first, hi=lo;
    for (const auto& [dm, rt] : rho)
    {
        lo.x=std::min(lo.x,dm.x); hi.x=std::max(hi.x,dm.x);
        lo.y=std::min(lo.y,dm.y); hi.y=std::max(hi.y,dm.y);
        lo.z=std::min(lo.z,dm.z); hi.z=std::max(hi.z,dm.z);
    }
    return {recip.GetCell().ToCartesian(rvec3_t(1,0,0)), recip.GetCell().ToCartesian(rvec3_t(0,1,0)),
            recip.GetCell().ToCartesian(rvec3_t(0,0,1)), lo, hi};
}

// Load once, then one value per point: each ro[p] is ONE thread's work, so the result is bit-identical at any thread count.
// It needs to be: under rho-tilde mixing the Fock build sees this density on every iteration.  Opt-in threading
// (QCHEM_OPENMP_THREADS; see qchem.Parallel).
static rvec_t SampleWith(FourierSampler& s, const ΔG_Map& rho, size_t n, double scale)
{
    s.Load(FourierTerms(rho.begin(), rho.end()));
    rvec_t ro(n);
    auto at=[&](size_t p) {ro[p]=scale*s.ValueAt(p);};
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::InnerThreads(); nthreads>1)
    {
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (size_t p=0; p<n; p++) at(p);
        return ro;
    }
#endif
    for (size_t p=0; p<n; p++) at(p);
    return ro;
}

rvec_t FourierMixCD::SampleOn(const PointSampling& ps) const
{
    const size_t n=ps.Points().size();
    if (itsRho.empty()) return rvec_t(n, 0.0);
    auto s=ps.MakeFourierSampler(BoxOf(itsRho, itsRecip));
    return SampleWith(*s, itsRho, n, itsScale);
}

rvec_t FourierMixCD::operator()(const rvec3vec_t& r) const
{
    if (itsRho.empty()) return rvec_t(r.size(), 0.0);
    auto s=MakeNUFFTSampler(r, BoxOf(itsRho, itsRecip), kNUFFTTolerance);
    return SampleWith(*s, itsRho, r.size(), itsScale);
}

// grad rho(r) = Sum_G i G rho-tilde(G) e^{iG.r} -- a two-line addition to the loop above, deliberately NOT
// written: nothing on the periodic path is a GGA or a gradient plotter, so an implementation here would be
// untested code.  Throwing keeps the FIRST such consumer honest; returning a silent zero (the shape this
// replaced) would hand it a plausible-looking wrong field instead (R1.4).
rvec3_t FourierMixCD::Gradient(const rvec3_t&) const
{
    throw std::logic_error("FourierMixCD::Gradient: grad(rho) is not implemented for a G-space density "
                           "(the periodic path is LDA-only).  Implement the i*G sum here, or ask the "
                           "density through a gradient-capable face.");
}

void FourierMixCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

} //namespace
