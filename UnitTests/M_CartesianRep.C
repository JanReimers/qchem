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
