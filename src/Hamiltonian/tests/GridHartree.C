// File: Hamiltonian/tests/GridHartree.C  The grid Hartree route: the direct Coulomb sum against a closed form.
//
// A Gaussian charge has a closed-form potential, q erf(sqrt(a) r)/r, so the sum over a mesh can be judged point by point
// without an SCF.  Tolerances are measured: the sum is a coarse rule, and these pin what it delivers, not an ideal.
#include "gtest/gtest.h"
#include <cmath>
#include <memory>
#include <numbers>
import qchem.Fitting.CoulombOperator;
import qchem.Blaze;
import qchem.Mesh;
import qchem.Structure;
import qchem.StructureData;
import qchem.Types;

using namespace qchem;

namespace {
const qcMesh::MeshParams theMesh = {.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=29, .beckeOrder=2};
qcMesh::Mesh OneAtomMesh()
{
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0,0,0)});
    return m.CreateIntegrationMesh(theMesh);
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
