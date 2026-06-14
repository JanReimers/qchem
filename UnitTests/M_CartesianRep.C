// File: UnitTests/M_CartesianRep.C  Cartesian-shell representation matrices (stage 2a).
#include <gtest/gtest.h>
#include <vector>
#include <array>
#include "blaze/Math.h"            // matrix operator* in this (non-module) TU
import qchem.Symmetry.CartesianRep;
import qchem.Symmetry.PointGroup;   // SymOp, to obtain operation matrices

using namespace Symmetry;

static std::vector<IVec3> Pshell() { return {{1,0,0},{0,1,0},{0,0,1}}; }
static std::vector<IVec3> Dshell() { return {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}}; }

// For a p-shell the representation matrix is the operation matrix R itself (p ~ x,y,z).
TEST(CartesianRep, p_shell_equals_R)
{
    rvec3_t z(0,0,1), y(0,1,0);
    for (const SymOp& op : { SymOp::Cn(z,2), SymOp::Cn(z,4), SymOp::Sigma(y), SymOp::Inversion() })
    {
        rmat_t D = CartesianShellRep(op.Matrix(), Pshell());
        for (int i=1;i<=3;i++) for (int j=1;j<=3;j++)
            EXPECT_NEAR(D(i-1,j-1), op.Matrix()(i,j), 1e-12);
    }
}

TEST(CartesianRep, identity_is_I)
{
    rmat_t D = CartesianShellRep(SymOp::E().Matrix(), Dshell());
    for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++)
        EXPECT_NEAR(D(i,j), (i==j)?1.0:0.0, 1e-12);
}

// Faithful representation: D(R1) D(R2) = D(R1 R2), checked on the d-shell.
TEST(CartesianRep, homomorphism_d_shell)
{
    rvec3_t z(0,0,1), x(1,0,0);
    Matrix3D<double> R1 = SymOp::Cn(z,4).Matrix();      // 90 deg about z
    Matrix3D<double> R2 = SymOp::Cn(x,2).Matrix();      // C2 about x
    rmat_t D1  = CartesianShellRep(R1,    Dshell());
    rmat_t D2  = CartesianShellRep(R2,    Dshell());
    rmat_t D12 = CartesianShellRep(R1*R2, Dshell());
    rmat_t prod = D1*D2;
    for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++)
        EXPECT_NEAR(prod(i,j), D12(i,j), 1e-10);
}

// d_{z2}-like component: under C2(z), every d monomial has even parity in (x,y) and is even
// in z, so the whole d-shell is invariant (D = I) for C2(z).
TEST(CartesianRep, d_shell_C2z_is_identity)
{
    rvec3_t z(0,0,1);
    rmat_t D = CartesianShellRep(SymOp::Cn(z,2).Matrix(), Dshell());
    // C2(z): x->-x, y->-y, z->z.  xx,yy,zz,xy unchanged; xz->-xz, yz->-yz.
    std::vector<double> diag = {1,1,1,1,-1,-1};   // matches Dshell() order
    for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++)
        EXPECT_NEAR(D(i,j), (i==j)?diag[i]:0.0, 1e-12);
}

// --- Full-basis representation on H2O (stage 2b) --------------------------------------
//
// AO layout (6 functions): O s [0], O p [1,2,3], H1 s [4], H2 s [5].  C2v with C2 along z.
static std::vector<AoShell> WaterShells()
{
    rvec3_t O(0,0,0.117), H1(0,0.757,-0.467), H2(0,-0.757,-0.467);
    return {
        {0, O,  {{0,0,0}},                  {1.0},         0},
        {1, O,  {{1,0,0},{0,1,0},{0,0,1}},  {1.0,1.0,1.0}, 1},
        {2, H1, {{0,0,0}},                  {1.0},         4},
        {2, H2, {{0,0,0}},                  {1.0},         5},
    };
}
static rvec3_t WaterOrigin()
{
    rvec3_t O(0,0,0.117), H1(0,0.757,-0.467), H2(0,-0.757,-0.467);
    return (O+H1+H2)/3.0;
}

TEST(CartesianRep, water_identity)
{
    rmat_t M = BuildOperationRep(WaterShells(), SymOp::E().Matrix(), WaterOrigin(), 1e-6);
    for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++)
        EXPECT_NEAR(M(i,j), (i==j)?1.0:0.0, 1e-12);
}

TEST(CartesianRep, water_C2_permutes_and_signs)
{
    rvec3_t z(0,0,1);
    rmat_t M = BuildOperationRep(WaterShells(), SymOp::Cn(z,2).Matrix(), WaterOrigin(), 1e-6);
    EXPECT_NEAR(M(0,0),  1.0, 1e-12);   // O s fixed
    EXPECT_NEAR(M(1,1), -1.0, 1e-12);   // O px -> -px
    EXPECT_NEAR(M(2,2), -1.0, 1e-12);   // O py -> -py
    EXPECT_NEAR(M(3,3),  1.0, 1e-12);   // O pz ->  pz
    EXPECT_NEAR(M(5,4),  1.0, 1e-12);   // H1 s -> H2 s
    EXPECT_NEAR(M(4,5),  1.0, 1e-12);   // H2 s -> H1 s
    EXPECT_NEAR(M(4,4),  0.0, 1e-12);   // (H s functions are swapped, not fixed)
}

// Representation law on the full H2O basis: M(g1) M(g2) = M(g1 g2).
TEST(CartesianRep, water_homomorphism)
{
    rvec3_t z(0,0,1), x(1,0,0);
    auto sh = WaterShells(); rvec3_t o = WaterOrigin();
    Matrix3D<double> R1 = SymOp::Cn(z,2).Matrix();    // C2(z)
    Matrix3D<double> R2 = SymOp::Sigma(x).Matrix();   // sigma(yz)
    rmat_t M1  = BuildOperationRep(sh, R1,    o, 1e-6);
    rmat_t M2  = BuildOperationRep(sh, R2,    o, 1e-6);
    rmat_t M12 = BuildOperationRep(sh, R1*R2, o, 1e-6);
    rmat_t prod = M1*M2;
    for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++)
        EXPECT_NEAR(prod(i,j), M12(i,j), 1e-10);
}
