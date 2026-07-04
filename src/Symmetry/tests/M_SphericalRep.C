// File: Symmetry/tests/M_SphericalRep.C  Real-spherical-harmonic operation reps (S1 of the spherical-SALC plan).
//
// SphericalShellRep projects the Cartesian rep through the harmonics' c2s expansion.  These tests pin the
// math independently of any basis: a p-shell in (x,y,z) order must reproduce R exactly, and the d-shell
// must be a faithful representation D(R1)D(R2)=D(R1R2) with the right C2(z) parities.
#include <gtest/gtest.h>
#include <vector>
#include <array>
#include <utility>
import qchem.Symmetry.Molecule.SphericalRep;
import qchem.Symmetry.Molecule.PointGroup;   // SymOp -> operation matrices
import qchem.Blaze;                  // rmat_t operator* in this (non-module) TU
using namespace qchem;
using namespace qchem::Symmetry;
using namespace qchem::Symmetry::Molecule;

// p harmonics ordered (x,y,z): each IS a Cartesian monomial, so the rep must equal R.
static HarmonicC2S Pxyz()
{
    return { {{{1,0,0},1.0}}, {{{0,1,0},1.0}}, {{{0,0,1},1.0}} };
}
// The five real d harmonics, m = -2..+2: xy, yz, (2z^2-x^2-y^2), xz, (x^2-y^2).
static HarmonicC2S Dshell()
{
    return {
        {{{1,1,0}, 1.0}},                                   // xy
        {{{0,1,1}, 1.0}},                                   // yz
        {{{0,0,2}, 2.0}, {{2,0,0},-1.0}, {{0,2,0},-1.0}},   // 2z^2-x^2-y^2
        {{{1,0,1}, 1.0}},                                   // xz
        {{{2,0,0}, 1.0}, {{0,2,0},-1.0}},                   // x^2-y^2
    };
}

// p-shell in (x,y,z) order: D_sph(R) == R for every operation.
TEST(SphericalRep, p_shell_equals_R)
{
    rvec3_t z(0,0,1), y(0,1,0), d(1,1,0);
    for (const SymOp& op : { SymOp::Cn(z,2), SymOp::Cn(z,4), SymOp::Cn(d,2), SymOp::Sigma(y), SymOp::Inversion() })
    {
        rmat_t D = SphericalShellRep(Pxyz()).Rep(op.Matrix());
        for (int i=1;i<=3;i++) for (int j=1;j<=3;j++)
            EXPECT_NEAR(D(i-1,j-1), op.Matrix()(i,j), 1e-12);
    }
}

TEST(SphericalRep, d_shell_identity_is_I)
{
    rmat_t D = SphericalShellRep(Dshell()).Rep(SymOp::E().Matrix());
    for (size_t i=0;i<5;i++) for (size_t j=0;j<5;j++)
        EXPECT_NEAR(D(i,j), (i==j)?1.0:0.0, 1e-12);
}

// Faithful representation on the d-shell: D(R1) D(R2) = D(R1 R2).
TEST(SphericalRep, homomorphism_d_shell)
{
    rvec3_t z(0,0,1), x(1,0,0);
    Matrix3D<double> R1 = SymOp::Cn(z,4).Matrix();      // 90 deg about z
    Matrix3D<double> R2 = SymOp::Cn(x,2).Matrix();      // C2 about x
    rmat_t D1  = SphericalShellRep(Dshell()).Rep(R1);
    rmat_t D2  = SphericalShellRep(Dshell()).Rep(R2);
    rmat_t D12 = SphericalShellRep(Dshell()).Rep(R1*R2);
    rmat_t prod = D1*D2;
    for (size_t i=0;i<5;i++) for (size_t j=0;j<5;j++)
        EXPECT_NEAR(prod(i,j), D12(i,j), 1e-10);
}

// Under C2(z) (x->-x, y->-y, z->z): xy(+) yz(-) z2(+) xz(-) x2-y2(+) -- a diagonal of those parities.
TEST(SphericalRep, d_shell_C2z_parities)
{
    rvec3_t z(0,0,1);
    rmat_t D = SphericalShellRep(Dshell()).Rep(SymOp::Cn(z,2).Matrix());
    std::vector<double> diag = {1,-1,1,-1,1};            // m = -2..+2 order of Dshell()
    for (size_t i=0;i<5;i++) for (size_t j=0;j<5;j++)
        EXPECT_NEAR(D(i,j), (i==j)?diag[i]:0.0, 1e-12);
}
