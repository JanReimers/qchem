// File: Hamiltonian/tests/GridHartree.C  The grid Hartree route: the direct Coulomb sum against a closed form.
//
// A Gaussian charge has a closed-form potential, q erf(sqrt(a) r)/r, so the sum over a mesh can be judged point by point
// without an SCF.  Tolerances are measured: the sum is a coarse rule, and these pin what it delivers, not an ideal.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
import qchem.Fitting.CoulombOperator;
import qchem.Blaze;
import qchem.Mesh;
import qchem.Mesh.Radial;
import qchem.Structure;
import qchem.StructureData;
import qchem.Types;

using namespace qchem;

namespace {
const qcMesh::MeshParams theMesh = {.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=29, .beckeOrder=2};
qcMesh::Mesh OneAtomMesh(const qcMesh::MeshParams& mp=theMesh)
{
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0,0,0)});
    return m.CreateIntegrationMesh(mp);
}
double Rho(double a, double r2) {return std::pow(a/std::numbers::pi,1.5)*std::exp(-a*r2);}   // unit charge
}

// The charge-weighted potential is the self-energy, which the erf form gives in closed form: E = q^2 sqrt(2a/pi).
TEST(DirectSumCoulomb, GaussianChargeSelfEnergy)
{
    const qcMesh::Mesh mesh=OneAtomMesh();
    const double a=1.3;
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const auto& r=mesh.Points()[g];
        q[g]=mesh.Weights()[g]*Rho(a,r.x*r.x+r.y*r.y+r.z*r.z);
    }
    const rvec_t v=Fitting::MakeDirectSumCoulomb(mesh)->Potential(q);
    double selfEnergy=0.0, charge=0.0;
    for (size_t g=0; g<mesh.size(); g++) {selfEnergy+=q[g]*v[g]; charge+=q[g];}
    EXPECT_NEAR(charge, 1.0, 1e-8);
    EXPECT_NEAR(selfEnergy, std::sqrt(2*a/std::numbers::pi), 1e-2);
}

// The radial Poisson solve on the same charge, against the closed form at every point: a one-site density is all l=0.  The error
// is the radial rule's, so it must fall with the radial count.
TEST(RadialPoissonCoulomb, GaussianChargePotential)
{
    auto worstAt=[](int nRadial)
    {
        qcMesh::MeshParams mp=theMesh; mp.nRadial=nRadial;
        const qcMesh::Mesh mesh=OneAtomMesh(mp);
        const double a=1.3;
        rvec_t q(mesh.size());
        for (size_t g=0; g<mesh.size(); g++)
        {
            const auto& r=mesh.Points()[g];
            q[g]=mesh.Weights()[g]*Rho(a,r.x*r.x+r.y*r.y+r.z*r.z);
        }
        const rvec_t v=Fitting::MakeRadialPoissonCoulomb(mesh, {rvec3_t(0,0,0)}, qcMesh::MakeRadial(mp).R(), 4)->Potential(q);
        double worst=0.0;
        for (size_t g=0; g<mesh.size(); g++)
        {
            const double r=norm(mesh.Points()[g]);
            if (r==0.0) continue;   // the rule's zero-weight node
            if (!std::isfinite(v[g])) return 1e300;
            worst=std::max(worst,std::fabs(v[g]-std::erf(std::sqrt(a)*r)/r));
        }
        return worst;
    };
    const double w80=worstAt(80), w160=worstAt(160), w320=worstAt(320);
    EXPECT_LT(w80,5e-5);
    EXPECT_LT(w160,w80/4);
    EXPECT_LT(w320,w160/4);
}

// A Gaussian off the centre needs every l: the expansion about the origin must converge on the closed form at the points.
TEST(RadialPoissonCoulomb, OffCentreGaussianNeedsHigherL)
{
    const qcMesh::Mesh mesh=OneAtomMesh();
    const double a=2.0; const rvec3_t c(0.3,-0.2,0.4);
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++)
    {
        const rvec3_t d=mesh.Points()[g]-c;
        q[g]=mesh.Weights()[g]*Rho(a,d*d);
    }
    auto worstAt=[&](int l)
    {
        const rvec_t v=Fitting::MakeRadialPoissonCoulomb(mesh, {rvec3_t(0,0,0)}, qcMesh::MakeRadial(theMesh).R(), l)->Potential(q);
        double worst=0.0;
        for (size_t g=0; g<mesh.size(); g+=13)
        {
            const double r=norm(mesh.Points()[g]-c);
            if (!std::isfinite(v[g])) return 1e300;
            worst=std::max(worst,std::fabs(v[g]-std::erf(std::sqrt(a)*r)/r));
        }
        return worst;
    };
    const double w2=worstAt(2), w10=worstAt(10);
    EXPECT_LT(w10,w2/20);
    EXPECT_LT(w10,1e-4);
}

// Linear in the charges, self term included.
TEST(DirectSumCoulomb, Linear)
{
    const qcMesh::Mesh mesh=OneAtomMesh();
    rvec_t q(mesh.size());
    for (size_t g=0; g<mesh.size(); g++) q[g]=mesh.Weights()[g]*std::exp(-double(g%7));
    const auto op=Fitting::MakeDirectSumCoulomb(mesh);
    const rvec_t v=op->Potential(q), v2=op->Potential(2.0*q);
    for (size_t g=0; g<mesh.size(); g+=997) EXPECT_NEAR(v2[g], 2.0*v[g], 1e-10*std::fabs(v[g]));
}
