// File: Math/Imp/PointSampling.C  The Fourier samplers a point set makes: the type-2 NUFFT (any point set) and the inverse FFT
// (the uniform midpoint raster).  The interface (Math/PointSampling.C) declares the faces; nothing outside names these classes.
module;
#include <complex>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <array>
#include <mutex>
#include <algorithm>
module qchem.PointSampling;
import qchem.NUFFT;
import qchem.FFT;

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
