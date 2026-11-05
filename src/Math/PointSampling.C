// File: Math/PointSampling.C  A point set that answers how a Fourier series is best evaluated on its own points.
//
// WHY (doc/OpenWork.md §4, the rho-tilde sampler; user rulings 2026-11-04): a density held as Fourier coefficients
// is sampled on a quadrature mesh every Kerker/Pulay iteration.  HOW is a property of the POINTS, not of the
// field and not of its client:
//   - any point set (Becke and every atom-centred mesh): a type-2 NUFFT (qchem.NUFFT), the series to a stated
//     tolerance;
//   - the uniform midpoint raster of the cell: a plain inverse FFT, EXACT and cheapest.
// So the point set MAKES the sampler (one virtual, NUFFT by default), and no caller asks what kind of mesh it holds.
// The field's half of the double dispatch is ScalarFunction<double>::SampleOn (qchem.ScalarFunction): a Fourier
// field overrides it to ask the point set for a sampler, and every other field evaluates itself at Points().
//
// A sampler is LOAD-ONCE, VALUE-PER-POINT: Load() does the per-coefficient-set work (deconvolve + FFT), ValueAt(p) is
// const and independent per point, so the caller threads the point loop (this leaf library has no threading) and each
// value is one thread's work -- bit-identical at any thread count.
module;
#include <complex>
#include <cmath>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <array>
#include <mutex>
#include <algorithm>
export module qchem.PointSampling;
import qchem.Types;   // rvec3_t, rvec3vec_t, ivec3_t, cvec_t, dcmplx
import qchem.NUFFT;
import qchem.FFT;
import qchem.CMath;   // Pi

export namespace qchem
{

//! The Fourier coefficients of a real-space field: f(r) = Sum_dm c(dm) e^{i G(dm).r}.
using FourierTerms = std::vector<std::pair<ivec3_t,dcmplx>>;
//! The reciprocal basis the series is written in, G(dm) = dm.x b1 + dm.y b2 + dm.z b3 (Cartesian, a.u.), and the
//! index box every dm of the series lies in.
struct FourierBox
{
    rvec3_t b1, b2, b3;
    ivec3_t lo, hi;
};

//! The series' REAL part on a fixed point set.  \c Load once per coefficient set, then \c ValueAt per point.
class FourierSampler
{
public:
    virtual ~FourierSampler() = default;
    virtual void   Load   (const FourierTerms&) =0;
    virtual double ValueAt(size_t p) const      =0;   //!< Re f(r_p), in the point set's own order
    virtual std::string Method() const          =0;   //!< what samples (the run report, and the gates that the right one answered)
};

//! \brief THE ACCURACY a point set's Fourier sampling is sized to (step 2, user rulings 2026-11-04): what its CONSUMER needs, stated by
//! the run (\c SolidCalculation) through the XC mesh recipe.  Default-constructed = UNSET: the NUFFT then runs at the fixed step-1
//! tolerance \c kNUFFTTolerance.
struct SamplingAccuracy
{
    double target         = 0.0;    //!< Ha: eps x atoms (an energy run) or eps (per matrix element, a geometry run); 0 = unset
    bool   perElement     = false;  //!< the run's client consumes matrix elements (a geometry), not only the energy
    double residualGate   = 0.0;    //!< the SCF's density gate (MinΔρPerElectron) the sampled residual pair is judged on; 0 = none
    double residualMargin = 0.1;    //!< kappa_R: the share of that gate the sampling error may use (GridAccuracy::samplerResidualMargin)
    double eta            = 1.0;    //!< the calibration delivered/bound the width choice divides by (GridAccuracy::samplerEta); 1 = rigorous
    //! The run's ONE tightening (<= 1): set by \c SolidCalculation at an SCF BOUNDARY after a previous SCF of the run measured the
    //! sampling over its target (the samp/τ column > 1), never inside an SCF.  1 = untightened.
    double tauScale       = 1.0;
};

//! \brief The point set's sampling SLOT: the accuracy it was built for, the ONE absolute target tau its consumer derives from it (fixed
//! from the seed density before the sampler's first sized use, and never changed inside the SCF), and what the sampler delivered.
//! Shared by every sampler the point set makes.  The sampler only MEASURES (an audit of every Load against the direct series on a point
//! subsample); it never changes its own settings in response -- an over-target is REPORTED (the samp/τ column) and acted on at the next
//! SCF boundary (user ruling 2026-11-04).
struct SamplingTarget
{
    SamplingAccuracy accuracy;
    double      tau     = 0.0;     //!< the SIZING tolerance: the width meets eta Bound ||c||_1 <= tau; 0 = unset (the step-1 fixed tolerance)
    //! The CONSUMER'S target max|delta rho| (before the run's one tightening, \c SamplingAccuracy::tauScale): what the delivered error is
    //! judged against.  The two differ only after a tightening -- which shrinks the sizing tau because the calibration proved optimistic,
    //! and cannot change what the consumer needs (judging against the tightened tau would read an over-optimistic eta as over target forever).
    double      tauTarget = 0.0;
    std::string tauFrom;           //!< which consumer bound set tau (the report)
    size_t      loads   = 0;       //!< sized Fourier samplings so far (tau is frozen once this is > 0)
    int         wMin    = 99, wMax = 0;
    double      worstDeliveredOverTau = 0.0;   //!< max over the slot's life of the audited delivered error / tauTarget
    double      windowWorst = 0.0;             //!< ...since the consumer last read it (\c TakeWindow): one SCF iteration's samp/τ
    //! The worst delivered/tau since the last call, and reset -- the consumer's per-iteration read.
    double TakeWindow() {const double w=windowWorst; windowWorst=0.0; return w;}
};

//! \brief A set of points that knows how a Fourier series is best evaluated on it.
class PointSampling
{
public:
    virtual ~PointSampling() = default;
    virtual const rvec3vec_t& Points() const =0;
    //! The sampling slot this point set's samplers size themselves by (null = none: the step-1 fixed tolerance).
    virtual SamplingTarget* FourierTarget() const {return nullptr;}
    //! The DEFAULT, right on ANY point set: a type-2 NUFFT at \c kNUFFTTolerance.  A point set with a structure the
    //! series can exploit (the uniform raster) overrides.
    virtual std::unique_ptr<FourierSampler> MakeFourierSampler(const FourierBox&) const;
};

//! The NUFFT's relative tolerance against ||c||_1.  STEP 1 (user, 2026-11-03): near machine precision, so replacing the
//! direct G-sum moves no anchor; sizing it from the run's eps (a rigorous bound + a measured calibration, CLAUDE.md
//! Design, pin 31) is step 2, its own increment.
inline constexpr double kNUFFTTolerance=1e-13;

//! The NUFFT sampler on \a points (any set) at the fixed relative \a tol (bare points: no consumer to size by).
std::unique_ptr<FourierSampler> MakeNUFFTSampler(const rvec3vec_t& points, const FourierBox&, double tol);
//! The NUFFT sampler on \a points sized by \a target (see \c SamplingTarget); a null or unset target is the fixed step-1 tolerance.
std::unique_ptr<FourierSampler> MakeNUFFTSampler(const rvec3vec_t& points, const FourierBox&, SamplingTarget* target);

//! \brief How a point set's builder records a structure the series can exploit -- the uniform override, HELD rather than
//! subclassed (a \c qcMesh::Mesh is a value type, returned and copied by value: a subclass would slice).
class FourierSamplingStrategy
{
public:
    virtual ~FourierSamplingStrategy() = default;
    virtual std::unique_ptr<FourierSampler> Make(const rvec3vec_t& points, const FourierBox&) const =0;
};

//! \brief The cell's uniform MIDPOINT raster, s = ((i+1/2)/n_1, (j+1/2)/n_2, (k+1/2)/n_3), points in (i, j, k) row-major
//! order (k fastest) -- what \c UnitCell::CreateIntegrationMesh builds.  The inverse FFT is EXACT for ANY n: the phase
//! e^{2 pi i dm.s} is n-periodic in dm, so each term lands in bin (dm mod n) carrying its half-voxel phase -- the same
//! series, no tolerance.  The series' b vectors must be the reciprocal of the raster's own cell.
class MidpointRasterSampling : public FourierSamplingStrategy
{
public:
    //! \a a1..a3 the raster's cell vectors: a series written on a DIFFERENT lattice (a supercell's) is refused, not mis-sampled.
    MidpointRasterSampling(const ivec3_t& n, const rvec3_t& a1, const rvec3_t& a2, const rvec3_t& a3) : itsN(n), itsA{a1,a2,a3} {}
    std::unique_ptr<FourierSampler> Make(const rvec3vec_t& points, const FourierBox&) const override;
    const ivec3_t& Divisions() const {return itsN;}
private:
    ivec3_t itsN;
    rvec3_t itsA[3];
};

} //namespace

namespace qchem
{
namespace
{
//! The plan's rigorous bound (NUFFT::Type2::Bound), cached: it depends only on the width and the box, which a run repeats.
double CachedBound(const ivec3_t& lo, const ivec3_t& hi, int w)
{
    static std::mutex m;
    static std::map<std::array<int,7>,double> cache;
    const std::array<int,7> key{w, lo.x,lo.y,lo.z, hi.x,hi.y,hi.z};
    {
        std::lock_guard<std::mutex> g(m);
        if (auto it=cache.find(key); it!=cache.end()) return it->second;
    }
    const double B=NUFFT::Type2::OfWidth(lo,hi,w).Bound();
    std::lock_guard<std::mutex> g(m);
    return cache[key]=B;
}

class NUFFTSampler : public FourierSampler
{
public:
    NUFFTSampler(const rvec3vec_t& r, const FourierBox& b, SamplingTarget* tgt, double tol)
        : itsBox(b), itsTarget(tgt), itsTol(tol), itsR(&r), itsS(r.size())
    {
        const double k=0.5/Pi;   // s_a = b_a.r / 2 pi: G.r = 2 pi dm.s
        for (size_t p=0; p<r.size(); p++) itsS[p]=rvec3_t(k*(b.b1*r[p]), k*(b.b2*r[p]), k*(b.b3*r[p]));
    }
    //! THE ONE PLACE the rigorous bound and the calibration meet (CLAUDE.md Design): w = the smallest width with
    //! eta * Bound(w) * ||c||_1 <= tau -- Bound the plan's rigorous alias bound (\c NUFFT::Type2::Bound), eta the stated calibration
    //! (\c GridAccuracy::samplerEta), tau the consumer's (fixed from the seed).  Then the Load is AUDITED against the direct series on
    //! ~200 of the points and the delivered/tau ratio recorded -- MEASURED, never acted on here.  Unset tau: the fixed step-1 tolerance.
    void Load(const FourierTerms& t) override
    {
        if (!itsTarget || !(itsTarget->tau>0.0))
        {
            itsPlan=std::make_unique<NUFFT::Type2>(itsBox.lo, itsBox.hi, itsTol);
            itsGrid=itsPlan->Grid(t);
            return;
        }
        SamplingTarget& T=*itsTarget;
        double l1=0.0;
        for (const auto& kc : t) l1+=std::abs(kc.second);
        const double eta=T.accuracy.eta;
        int w=2;
        while (w<16 && eta*CachedBound(itsBox.lo,itsBox.hi,w)*l1>T.tau) w++;
        itsPlan=std::make_unique<NUFFT::Type2>(NUFFT::Type2::OfWidth(itsBox.lo, itsBox.hi, w));
        itsGrid=itsPlan->Grid(t);
        itsDelivered=Audit(t);
        const double r=itsDelivered/(T.tauTarget>0.0 ? T.tauTarget : T.tau);   // judged against what the consumer needs
        T.loads++; T.wMin=std::min(T.wMin,w); T.wMax=std::max(T.wMax,w);
        T.worstDeliveredOverTau=std::max(T.worstDeliveredOverTau, r);
        T.windowWorst=std::max(T.windowWorst, r);
    }
    double ValueAt(size_t p) const override {return std::real(itsPlan->Interpolate(itsGrid, itsS[p]));}
    std::string Method() const override
    {
        if (!itsPlan) return "NUFFT (not loaded)";   // the width is chosen per Load
        std::string m="NUFFT (w "+std::to_string(itsPlan->Width());
        if (itsTarget && itsTarget->tau>0.0) m+=", tau "+std::to_string(itsTarget->tau)+", delivered "+std::to_string(itsDelivered);
        return m+")";
    }
private:
    //! max |NUFFT - direct series| (real parts) on ~200 of the points: the direct sum by per-axis phase tables, one complex
    //! multiply per term per point.
    double Audit(const FourierTerms& t) const
    {
        const size_t n=itsS.size(), stride=std::max<size_t>(1, n/200);
        const ivec3_t lo=itsBox.lo, hi=itsBox.hi;
        std::vector<dcmplx> tx(hi.x-lo.x+1), ty(hi.y-lo.y+1), tz(hi.z-lo.z+1);
        double e=0.0;
        for (size_t p=0; p<n; p+=stride)
        {
            const rvec3_t& s=itsS[p];
            for (int m=lo.x; m<=hi.x; m++) tx[m-lo.x]=std::exp(dcmplx(0.0,2.0*Pi*m*s.x));
            for (int m=lo.y; m<=hi.y; m++) ty[m-lo.y]=std::exp(dcmplx(0.0,2.0*Pi*m*s.y));
            for (int m=lo.z; m<=hi.z; m++) tz[m-lo.z]=std::exp(dcmplx(0.0,2.0*Pi*m*s.z));
            dcmplx f(0.0);
            for (const auto& [dm, c] : t) f+=c*tx[dm.x-lo.x]*ty[dm.y-lo.y]*tz[dm.z-lo.z];
            e=std::max(e, std::abs(std::real(itsPlan->Interpolate(itsGrid, s))-std::real(f)));
        }
        return e;
    }
    FourierBox                     itsBox;
    SamplingTarget*                itsTarget;
    double                         itsTol;
    const rvec3vec_t*              itsR;
    std::vector<rvec3_t>           itsS;      //!< the points in fractional coordinates
    std::unique_ptr<NUFFT::Type2>  itsPlan;
    cvec_t                         itsGrid;
    double                         itsDelivered=0.0;
};

class RasterSampler : public FourierSampler
{
public:
    RasterSampler(const ivec3_t& n, size_t npts) : itsN(n)
    {
        if (npts!=size_t(n.x)*n.y*n.z)
            throw std::logic_error("MidpointRasterSampling: the point set is not the n1 x n2 x n3 raster it was installed for");
    }
    void Load(const FourierTerms& t) override
    {
        itsGrid=cvec_t(size_t(itsN.x)*itsN.y*itsN.z, dcmplx(0.0));
        const size_t n1=itsN.y, n2=itsN.z;
        for (const auto& [dm, c] : t)
        {
            const double half=Pi*(double(dm.x)/itsN.x+double(dm.y)/itsN.y+double(dm.z)/itsN.z);   // the half-voxel offset
            itsGrid[(Wrap(dm.x,itsN.x)*n1+Wrap(dm.y,itsN.y))*n2+Wrap(dm.z,itsN.z)] += c*std::exp(dcmplx(0.0,half));
        }
        itsGrid=FFT::FFT3D(std::move(itsGrid), itsN, +1);
    }
    double ValueAt(size_t p) const override {return std::real(itsGrid[p]);}
    std::string Method() const override {return "inverse FFT (midpoint raster)";}
private:
    static size_t Wrap(int j, int n) {const int m=j%n; return size_t(m<0 ? m+n : m);}
    ivec3_t itsN;
    cvec_t  itsGrid;
};
} //anonymous

std::unique_ptr<FourierSampler> MakeNUFFTSampler(const rvec3vec_t& points, const FourierBox& b, double tol)
{
    return std::make_unique<NUFFTSampler>(points, b, nullptr, tol);
}
std::unique_ptr<FourierSampler> MakeNUFFTSampler(const rvec3vec_t& points, const FourierBox& b, SamplingTarget* target)
{
    return std::make_unique<NUFFTSampler>(points, b, target, kNUFFTTolerance);
}

std::unique_ptr<FourierSampler> PointSampling::MakeFourierSampler(const FourierBox& b) const
{
    return MakeNUFFTSampler(Points(), b, FourierTarget());
}

std::unique_ptr<FourierSampler> MidpointRasterSampling::Make(const rvec3vec_t& points, const FourierBox& b) const
{
    const rvec3_t* bs[3]={&b.b1,&b.b2,&b.b3};
    for (int i=0; i<3; i++) for (int j=0; j<3; j++)
        if (std::abs(itsA[i]*(*bs[j])-(i==j ? 2.0*Pi : 0.0))>1e-9*2.0*Pi)
            throw std::logic_error("MidpointRasterSampling: the series is not written on this raster's lattice "
                                   "(a_i.b_j != 2 pi delta_ij) -- its inverse FFT would sample a different function");
    return std::make_unique<RasterSampler>(itsN, points.size());
}

} //namespace
