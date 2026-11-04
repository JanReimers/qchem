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

//! \brief A set of points that knows how a Fourier series is best evaluated on it.
class PointSampling
{
public:
    virtual ~PointSampling() = default;
    virtual const rvec3vec_t& Points() const =0;
    //! The DEFAULT, right on ANY point set: a type-2 NUFFT at \c kNUFFTTolerance.  A point set with a structure the
    //! series can exploit (the uniform raster) overrides.
    virtual std::unique_ptr<FourierSampler> MakeFourierSampler(const FourierBox&) const;
};

//! The NUFFT's relative tolerance against ||c||_1.  STEP 1 (user, 2026-11-03): near machine precision, so replacing the
//! direct G-sum moves no anchor; sizing it from the run's eps (a rigorous bound + a measured calibration, CLAUDE.md
//! Design, pin 31) is step 2, its own increment.
inline constexpr double kNUFFTTolerance=1e-13;

//! The NUFFT sampler on \a points (any set).
std::unique_ptr<FourierSampler> MakeNUFFTSampler(const rvec3vec_t& points, const FourierBox&, double tol);

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
class NUFFTSampler : public FourierSampler
{
public:
    NUFFTSampler(const rvec3vec_t& r, const FourierBox& b, double tol) : itsPlan(b.lo, b.hi, tol), itsS(r.size())
    {
        const double k=0.5/Pi;   // s_a = b_a.r / 2 pi: G.r = 2 pi dm.s
        for (size_t p=0; p<r.size(); p++) itsS[p]=rvec3_t(k*(b.b1*r[p]), k*(b.b2*r[p]), k*(b.b3*r[p]));
    }
    void   Load   (const FourierTerms& t) override {itsGrid=itsPlan.Grid(t);}
    double ValueAt(size_t p) const override {return std::real(itsPlan.Interpolate(itsGrid, itsS[p]));}
    std::string Method() const override {return "NUFFT (w "+std::to_string(itsPlan.Width())+")";}
private:
    NUFFT::Type2         itsPlan;
    std::vector<rvec3_t> itsS;      //!< the points in fractional coordinates
    cvec_t               itsGrid;
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
    return std::make_unique<NUFFTSampler>(points, b, tol);
}

std::unique_ptr<FourierSampler> PointSampling::MakeFourierSampler(const FourierBox& b) const
{
    return MakeNUFFTSampler(Points(), b, kNUFFTTolerance);
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
