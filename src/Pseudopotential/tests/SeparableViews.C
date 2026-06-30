// File: src/Pseudopotential/tests/SeparableViews.C  The two radial views of a separable PP are consistent.
//
// SeparablePotential exposes a projector's radial function in two spectral views: the reciprocal leaf
// Projector(q) (plane-wave) and the real leaf BetaR(r) (molecular/atomic).  They MUST be the same function,
// related by the spherical-Bessel transform:
//     integral_0^inf BetaR_p(r) j_l(q r) r^2 dr  ==  Projector_p(q) / sqrt(4 pi).
// The 1/sqrt(4pi) is fixed by KB consistency: <plane wave|beta Y_lm> = (4pi/sqrt Omega) i^l Y*_lm(q) beta~(q),
// and summing |beta Y_lm><beta Y_lm| over m (addition theorem -> (2l+1)/4pi P_l) reproduces the reciprocal
// assembler's (2l+1) P_l Projector*Projector iff Projector = sqrt(4pi) beta~.  This nails the real-space
// normalisation (the KB energy is unforgiving) before any assembly relies on it.
#include "gtest/gtest.h"
#include <cmath>
#include <vector>
#include <initializer_list>

import qchem.Pseudopotential.GTH_Potentials;     // GetGTH, GTH_PP (real HGH parameters)
import qchem.Mesh.Quadrature;                     // qcMesh::RadialMesh, MakeRadial
import qchem.Math;                                // Pi
using namespace qchem;

using namespace qchem::Pseudopotential;

// integral_0^inf BetaR_p(r) j_l(q r) r^2 dr, on a fine log radial mesh (the weights fold in r^2).
static double BesselTransform(const HGH_SeparablePotential& sep, int Z, size_t p, int l, double q,
                              const qcMesh::RadialMesh& mesh)
{
    const rvec_t& R=mesh.R();  const rvec_t& W=mesh.W();
    double s=0;
    for (size_t i=0;i<R.size();i++)
        s += W[i] * sep.BetaR(Z,p,R[i]) * std::sph_bessel((unsigned)l, q*R[i]);
    return s;
}

TEST(SeparablePotentialViews, RealBesselTransformMatchesReciprocal)
{
    const int Z=14;                                  // Si
    GTH_PP pp = GetGTH("Si", "LDA", 4);              // real HGH q4 parameters (local + KB-separable)
    const HGH_SeparablePotential& sep = pp.nonlocal;
    ASSERT_GT(sep.NumProjectors(Z), 0u);

    // Fine log mesh: the HGH projectors are Gaussians (decay within a few r_loc ~ a couple of bohr).
    qcMesh::RadialMesh mesh = qcMesh::MakeRadial({.radial=qcMesh::RadialKind::Log, .nRadial=4000,
                                                  .logStart=1e-4, .logStop=12.0});
    const double s4pi = std::sqrt(4.0*Pi);

    for (size_t p=0;p<sep.NumProjectors(Z);++p)
    {
        int l = sep.AngularMomentum(Z,p);
        for (double q : {0.5, 1.0, 2.0, 3.0, 5.0})
        {
            double lhs = BesselTransform(sep, Z, p, l, q, mesh);   // real view, transformed: beta~(q)
            double rhs = sep.Projector(Z, p, q) / s4pi;            // reciprocal view: Projector/sqrt(4pi)
            EXPECT_NEAR(lhs, rhs, 1e-3*std::fabs(rhs) + 1e-5)
                << "projector p=" << p << " (l=" << l << ") at q=" << q;
        }
    }
}

// Real GTH parameters that actually carry an f channel: Cs q9 has l=0,1,2 AND an l=3 projector (as do the
// lanthanides/actinides).  Before the table grew past l=2 this channel hit assert(false); here the real
// f-projector's two radial views are checked consistent, end to end on database parameters.
TEST(SeparablePotentialViews, CesiumFChannelBesselTransformMatchesReciprocal)
{
    const int Z=55;
    GTH_PP pp = GetGTH("Cs", "LDA", 9);
    const HGH_SeparablePotential& sep = pp.nonlocal;
    int maxl=0; for (size_t p=0;p<sep.NumProjectors(Z);++p) maxl=std::max(maxl, sep.AngularMomentum(Z,p));
    ASSERT_EQ(maxl, 3) << "Cs q9 should carry an l=3 (f) projector";

    qcMesh::RadialMesh mesh = qcMesh::MakeRadial({.radial=qcMesh::RadialKind::Log, .nRadial=4000,
                                                  .logStart=1e-4, .logStop=12.0});
    const double s4pi = std::sqrt(4.0*Pi);
    for (size_t p=0;p<sep.NumProjectors(Z);++p)
    {
        int l = sep.AngularMomentum(Z,p);
        for (double q : {0.5, 1.0, 2.0, 3.0, 5.0})
            EXPECT_NEAR(BesselTransform(sep, Z, p, l, q, mesh), sep.Projector(Z, p, q)/s4pi,
                        1e-3*std::fabs(sep.Projector(Z,p,q)/s4pi) + 1e-5)
                << "Cs projector p=" << p << " (l=" << l << ") at q=" << q;
    }
}

// Exercise EVERY tabulated momentum-space polynomial Q_i^l -- in particular the higher projectors (l=1 i=2,
// l=2 i=1) and the f channel (l=3 i=0,1) that the analytic HGH table grew to cover.  A DIAGONAL h-matrix
// makes each KB projector a PURE single radial p_i^l (the eigenvectors are unit vectors), so the same
// real<->reciprocal Bessel-transform identity validates each Q_i^l individually.  (The real-space p_i^l is
// general in (l,i); this pins the matching closed-form Q_i^l, the part that was hand-tabulated.)
TEST(SeparablePotentialViews, AllChannelsBesselTransformMatchesReciprocal)
{
    const int Z=1;
    const double rl=0.42;
    auto diag=[](std::initializer_list<double> ds){
        std::vector<std::vector<double>> h(ds.size(), std::vector<double>(ds.size(), 0.0));
        size_t k=0; for (double d:ds) { h[k][k]=d; ++k; } return h;
    };
    HGH_SeparablePotential sep;
    sep.AddChannel(1, rl, diag({ 1.3, -0.7,  0.4}));   // l=1: i=0,1,2
    sep.AddChannel(2, rl, diag({-0.9,  0.5}));         // l=2: i=0,1
    sep.AddChannel(3, rl, diag({ 0.8, -0.6}));         // l=3: i=0,1 (the f channel)

    qcMesh::RadialMesh mesh = qcMesh::MakeRadial({.radial=qcMesh::RadialKind::Log, .nRadial=4000,
                                                  .logStart=1e-4, .logStop=12.0});
    const double s4pi = std::sqrt(4.0*Pi);
    ASSERT_EQ(sep.NumProjectors(Z), 7u);               // 3 + 2 + 2
    for (size_t p=0;p<sep.NumProjectors(Z);++p)
    {
        int l = sep.AngularMomentum(Z,p);
        for (double q : {0.5, 1.0, 2.0, 3.0, 5.0})
        {
            double lhs = BesselTransform(sep, Z, p, l, q, mesh);
            double rhs = sep.Projector(Z, p, q) / s4pi;
            EXPECT_NEAR(lhs, rhs, 1e-3*std::fabs(rhs) + 1e-5)
                << "synthetic projector p=" << p << " (l=" << l << ") at q=" << q;
        }
    }
}
