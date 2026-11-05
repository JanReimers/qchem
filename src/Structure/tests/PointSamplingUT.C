// File: src/Structure/tests/PointSamplingUT.C  The point set chooses how a Fourier series is sampled on it (qchem.PointSampling).
//
// The rho-tilde sampler (doc/OpenWork.md §4, user rulings 2026-11-04): a Becke mesh takes the NUFFT default, the cell's uniform
// midpoint raster overrides with the exact inverse FFT, and a field that is NOT a Fourier series is sampled exactly as before.
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
import qchem.Lattice_3D;      // UnitCell / FCCUnitCell
import qchem.Mesh;            // qcMesh::Mesh (a PointSampling)
import qchem.Mesh.XCPolicy;   // BeckeXCParams
import qchem.PointSampling;
import qchem.ScalarFunction;
import qchem.Types;
using namespace qchem;

namespace
{
UnitCell SiPrimitive()
{
    FCCUnitCell c(10.26);
    c.AddAtom(14, {0.00,0.00,0.00});
    c.AddAtom(14, {0.25,0.25,0.25});
    return c;
}
//! A ball of terms on \a cell's reciprocal lattice with NO conjugate symmetry, and its box and l1 norm.
struct Series {FourierTerms t; FourierBox box; double l1=0.0;};
Series Ball(const UnitCell& cell, int m)
{
    const CellGeometry rc=cell.MakeReciprocal();
    Series s;
    s.box={rc.ToCartesian(rvec3_t(1,0,0)), rc.ToCartesian(rvec3_t(0,1,0)), rc.ToCartesian(rvec3_t(0,0,1)), ivec3_t(-m,-m,-m), ivec3_t(m,m,m)};
    for (int x=-m; x<=m; x++) for (int y=-m; y<=m; y++) for (int z=-m; z<=m; z++)
        if (x*x+y*y+z*z<=m*m)
        {
            const dcmplx c(std::sin(1.3*x+0.7*y-0.4*z)/(1.0+x*x+y*y+z*z), std::cos(0.9*x-1.1*y+0.3*z)/(2.0+x*x+y*y+z*z));
            s.t.push_back({ivec3_t(x,y,z),c}); s.l1+=std::abs(c);
        }
    return s;
}
double Direct(const Series& s, const rvec3_t& r)
{
    dcmplx f(0.0);
    for (const auto& [dm,c] : s.t)
    {
        const rvec3_t G=double(dm.x)*s.box.b1+double(dm.y)*s.box.b2+double(dm.z)*s.box.b3;
        f+=c*std::exp(dcmplx(0.0, G*r));
    }
    return std::real(f);
}
//! The max deviation of \a mesh's own sampler from the direct series, and which sampler answered.
std::pair<double,std::string> Deviation(const qcMesh::Mesh& mesh, const Series& s)
{
    auto smp=mesh.MakeFourierSampler(s.box);
    smp->Load(s.t);
    double e=0.0;
    for (size_t p=0; p<mesh.size(); p+=7) e=std::max(e, std::abs(smp->ValueAt(p)-Direct(s,mesh.Points()[p])));
    return {e, smp->Method()};
}
struct Plain : ScalarFunction<double>   //!< a field that is not a Fourier series
{
    using ScalarFunction<double>::operator();   // keep the batched overload visible
    double operator()(const rvec3_t& r) const override {return std::sin(r.x)*std::cos(0.3*r.y)+r.z*r.z;}
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}
};
}

// A Becke (atom-centred) mesh takes the NUFFT default, to the series at its stated tolerance.
TEST(PointSampling, ABeckeMeshSamplesTheSeriesByNUFFT)
{
    const UnitCell cell=SiPrimitive();
    const qcMesh::Mesh m=cell.CreateIntegrationMesh(qcMesh::BeckeXCParams(10, 2.0, 11));
    const Series s=Ball(cell, 6);
    const auto [e, how]=Deviation(m, s);
    EXPECT_EQ(how.rfind("NUFFT",0), 0u) << how;
    EXPECT_LT(e, kNUFFTTolerance*s.l1);
}

// The uniform midpoint raster answers with the EXACT inverse FFT -- including a series wider than the raster's Nyquist,
// which the bin-folding samples exactly (the phase is n-periodic in dm) -- and a copy keeps that sampling.
TEST(PointSampling, TheUniformRasterSamplesTheSeriesByAnExactInverseFFT)
{
    const UnitCell cell=SiPrimitive();
    const qcMesh::Mesh m=cell.CreateIntegrationMesh({.nUniform=9, .cellKind=qcMesh::UnitCellKind::Uniform});
    for (int ball : {3, 7})   // 7: |dm| up to 7 on a 9-point raster
    {
        const Series s=Ball(cell, ball);
        const auto [e, how]=Deviation(m, s);
        EXPECT_EQ(how, "inverse FFT (midpoint raster)");
        EXPECT_LT(e, 1e-14*s.l1) << "ball " << ball;
    }
    const qcMesh::Mesh copy=m;
    EXPECT_EQ(copy.MakeFourierSampler(Ball(cell,3).box)->Method(), "inverse FFT (midpoint raster)");
}

// The raster strategy is a statement about the points: moving them drops it (the NUFFT default answers), and a series
// written on ANOTHER lattice is refused rather than sampled as a different function.
TEST(PointSampling, TheRasterStrategyGuardsItsOwnPointsAndLattice)
{
    const UnitCell cell=SiPrimitive();
    qcMesh::Mesh m=cell.CreateIntegrationMesh({.nUniform=6, .cellKind=qcMesh::UnitCellKind::Uniform});
    Series other=Ball(cell, 2);
    other.box.b1=0.5*other.box.b1;   // a supercell's reciprocal vector
    EXPECT_THROW(m.MakeFourierSampler(other.box), std::logic_error);
    m.ShiftOrigin(rvec3_t(0.1,0,0));
    EXPECT_EQ(m.MakeFourierSampler(Ball(cell,2).box)->Method().rfind("NUFFT",0), 0u);
}

// A field that is NOT a Fourier series is sampled on any point set exactly as f(points) -- BITWISE, on both meshes.
TEST(PointSampling, ANonFourierFieldIsSampledBitwiseAsBefore)
{
    const UnitCell cell=SiPrimitive();
    const Plain f;
    for (const qcMesh::Mesh& m : {cell.CreateIntegrationMesh({.nUniform=6, .cellKind=qcMesh::UnitCellKind::Uniform}),
                                  cell.CreateIntegrationMesh(qcMesh::BeckeXCParams(6, 2.0, 7))})
    {
        const rvec_t a=f.SampleOn(m), b=f(m.Points());
        ASSERT_EQ(a.size(), b.size());
        for (size_t p=0; p<a.size(); p++) ASSERT_EQ(a[p], b[p]) << "point " << p;
    }
}
