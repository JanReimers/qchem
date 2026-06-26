// File: src/Structure/tests/MolecularMeshTests.C  Acceptance tests for the Becke molecular mesh.
//
// The mirror of the old MoleculeMesh bug: a homonuclear dimer at separation -> 0 must integrate to
// the SINGLE-atom result, and a well-separated dimer must be ADDITIVE (= 2 x the atom).
#include "gtest/gtest.h"
#include <cmath>

import qchem.Structure;                 // Molecule, Atom, Structure
import qchem.Structure.MolecularMesh;  // MakeMolecularMesh
import qchem.Mesh.Quadrature;          // qcMesh::Integrate, ScalarField
import qchem.Math;                      // Pi32

// Note: this TU sees BOTH the old global `Mesh`/`MeshParams` (via qchem.Structure) and the new
// qcMesh::* -- so we qualify qcMesh:: throughout and never name a bare Mesh/MeshParams.

namespace
{
// Gaussian exp(-|r-c|^2); integral over R^3 = pi^{3/2}.
class GaussAt : public qcMesh::ScalarField<double>
{
    rvec3_t itsC;
public:
    explicit GaussAt(const rvec3_t& c) : itsC(c) {}
    double  operator()(const rvec3_t& r) const override {double m=norm(r-itsC); return std::exp(-m*m);}
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};

// Sum of two unit Gaussians centred at a and b; integral = 2 pi^{3/2}.
class TwoGauss : public qcMesh::ScalarField<double>
{
    rvec3_t itsA, itsB;
public:
    TwoGauss(const rvec3_t& a, const rvec3_t& b) : itsA(a), itsB(b) {}
    double  operator()(const rvec3_t& r) const override
    { double ma=norm(r-itsA), mb=norm(r-itsB); return std::exp(-ma*ma)+std::exp(-mb*mb); }
    rvec3_t Gradient  (const rvec3_t&)   const override {return rvec3_t(0,0,0);}
};

qcMesh::MeshParams Params()
{
    qcMesh::MeshParams mp;
    mp.radial=qcMesh::RadialKind::MHL; mp.nRadial=40; mp.mhl_m=2; mp.mhl_alpha=2.0;
    mp.angular=qcMesh::AngularKind::Gauss; mp.nAngular=24;
    return mp;
}
} //anon

// The bug's mirror: two coincident atoms must integrate to exactly the single-atom result.
TEST(MolecularMesh, CoincidentDimerEqualsAtom)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t o(0,0,0);

    Atom a(1,o);
    auto m1=MakeMolecularMesh(a,mp);

    Molecule dimer;
    dimer.Insert(new Atom(1,o));
    dimer.Insert(new Atom(1,o));          // right on top of each other (R_ab = 0)
    auto m2=MakeMolecularMesh(dimer,mp);

    GaussAt g(o);
    double I1=qcMesh::Integrate(m1,g);
    double I2=qcMesh::Integrate(m2,g);
    EXPECT_NEAR(I2,I1,1e-12);             // coincident dimer == single atom (no 0/0 corruption)
    EXPECT_NEAR(I1,Pi32,1e-4);            // and both reproduce integral exp(-r^2) d^3r = pi^{3/2}
}

// Additivity: a well-separated dimer integrates the sum of two atom-centred Gaussians to 2 pi^{3/2}.
TEST(MolecularMesh, SeparatedDimerIsAdditive)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t RA(-8,0,0), RB(8,0,0);

    Molecule dimer;
    dimer.Insert(new Atom(1,RA));
    dimer.Insert(new Atom(1,RB));
    auto m=MakeMolecularMesh(dimer,mp);

    TwoGauss f(RA,RB);
    EXPECT_NEAR(qcMesh::Integrate(m,f),2*Pi32,1e-3);
}

// A far second atom must not corrupt the integral of a function localised on the first.
TEST(MolecularMesh, FarAtomDoesNotCorrupt)
{
    qcMesh::MeshParams mp=Params();
    rvec3_t RA(-8,0,0), RB(8,0,0);

    Molecule dimer;
    dimer.Insert(new Atom(1,RA));
    dimer.Insert(new Atom(1,RB));
    auto m=MakeMolecularMesh(dimer,mp);

    GaussAt g(RA);                        // localised on atom A only
    EXPECT_NEAR(qcMesh::Integrate(m,g),Pi32,1e-3);
}
