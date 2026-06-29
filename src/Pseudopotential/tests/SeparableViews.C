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

import qchem.Pseudopotential.GTH_Potentials;     // GetGTH, GTH_PP (real HGH parameters)
import qchem.Mesh.Quadrature;                     // qcMesh::RadialMesh, MakeRadial
import qchem.Math;                                // Pi

using namespace Pseudopotential;

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
