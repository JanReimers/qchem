// File: Mesh/tests/IntegratorUT.C  MatrixIntegrator -- the forward/adjoint pair, and its adjointness.
//
// The class exists so a forward and an adjoint cannot be mismatched (doc/CleanupCandidates.md R1.0j).
// The property that makes the pairing WORTH enforcing is testable directly and is the main test here:
//
//     <v, Forward(D)>_w  ==  <Adjoint(v), D>        for every D and every v
//
// i.e. the two maps really are adjoint under the mesh's own inner product.  If a future realization
// screens one direction and not the other, THIS is what breaks -- which is exactly the defect that cost
// the Si run 14 -> 60 iterations before the two halves were put on one object.
#include "gtest/gtest.h"
#include <cmath>
#include <memory>

import qchem.Mesh.Integrator;
import qchem.Mesh.Builder;
import qchem.VectorFunction;
import qchem.Types;

using namespace qchem;

namespace {

//! Three fixed, linearly independent real functions -- enough to make D and v generic.  Geometry is
//! irrelevant: the adjoint identity is algebra over whatever the basis returns.
class TinyBasis : public virtual VectorFunction<double>
{
public:
    virtual size_t   GetVectorSize() const override {return 3;}
    virtual rvec_t   operator()(const rvec3_t& r) const override
    {
        const double x=r.x;
        return rvec_t{1.0, x, std::exp(-0.5*x*x)};
    }
    virtual rvec3vec_t Gradient(const rvec3_t&) const override {return rvec3vec_t(3);}
};

qcMesh::Mesh SpreadMesh(size_t n)
{
    qcMesh::MeshBuilder b;
    for (size_t i=0;i<n;i++) b.Append(rvec3_t(-2.0+4.0*double(i)/double(n-1),0,0), 0.3+0.1*double(i%3));
    return b.take();
}

} // anon

//---------------------------------------------------------------------------------------------------
// ★ THE PROPERTY THE CLASS EXISTS TO PROTECT.  <v, Forward(D)>_w == <Adjoint(v), D>, to machine
// precision, for a generic D and a generic v.  This is what "H = dE/dD" rests on.
TEST(MatrixIntegrator, ForwardAndAdjointAreExactlyAdjoint)
{
    const qcMesh::Mesh m=SpreadMesh(24);
    const TinyBasis    b;
    const qcMesh::DenseMatrixIntegrator<double> I(m, b);

    hmat_t<double> D(3);
    D(0,0)=1.5; D(0,1)=-0.4; D(0,2)=0.7;
            D(1,1)=2.0; D(1,2)=-1.1;
                    D(2,2)=0.9;

    rvec_t v(m.size());
    for (size_t g=0; g<m.size(); g++) v[g]=std::sin(1.7*double(g))+0.5;

    // LHS: integrate v against the collocated density, using the mesh's own quadrature
    const rvec_t rho=I.Forward(D);
    ASSERT_EQ(rho.size(), m.size());
    rvec_t vrho(m.size());
    for (size_t g=0; g<m.size(); g++) vrho[g]=v[g]*rho[g];
    const double lhs=I.Integrate(vrho);

    // RHS: contract the assembled matrix with D.  Full double sum, since D and the matrix are Hermitian
    // adaptors and the off-diagonals must count twice.
    const hmat_t<double> H=I.Adjoint(v);
    double rhs=0.0;
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) rhs += H(i,j)*D(i,j);

    EXPECT_NEAR(lhs, rhs, 1e-12*std::max(1.0,std::abs(lhs)))
        << "Forward and Adjoint are not adjoint -- the pairing this class exists to guarantee is broken";
}

//---------------------------------------------------------------------------------------------------
// A density matrix produces a REAL, and for a positive-semidefinite D a NON-NEGATIVE, density.  The XC
// path relies on rho>=0 pointwise (the rho>0 guard is inert on an aufbau D), so it is worth pinning.
TEST(MatrixIntegrator, ForwardOfAPSDDensityMatrixIsNonNegative)
{
    const qcMesh::Mesh m=SpreadMesh(16);
    const TinyBasis    b;
    const qcMesh::DenseMatrixIntegrator<double> I(m, b);

    hmat_t<double> D(3);                  // an outer product c c^T: PSD by construction
    const rvec_t c{0.6,-1.2,0.8};
    for (size_t i=0;i<3;i++) for (size_t j=i;j<3;j++) D(i,j)=c[i]*c[j];

    const rvec_t rho=I.Forward(D);
    for (size_t g=0; g<rho.size(); g++)
        EXPECT_GE(rho[g], -1e-14) << "a PSD density matrix must collocate to a non-negative density";
}

//---------------------------------------------------------------------------------------------------
// Integrate is the mesh's own rule, and NumPoints sizes both arrays -- so a caller can allocate a field
// without reaching past the interface for the mesh.
TEST(MatrixIntegrator, IntegrateIsTheMeshRuleAndNumPointsSizesTheArrays)
{
    const qcMesh::Mesh m=SpreadMesh(10);
    const TinyBasis    b;
    const qcMesh::DenseMatrixIntegrator<double> I(m, b);

    EXPECT_EQ(I.NumPoints(), m.size());
    rvec_t ones(I.NumPoints(), 1.0);
    double w=0.0; for (size_t g=0;g<m.size();g++) w+=m.Weights()[g];
    EXPECT_DOUBLE_EQ(I.Integrate(ones), w);      // integral of 1 == total weight
}

//---------------------------------------------------------------------------------------------------
// ★ THE ISP SPLIT (user, 2026-09-09): "we need to SOLID::ISP MatrixIntegrator so the two sides are split,
// but only at the abstract interface level.  Charge density holds pointer to the forward interface side,
// and basis holds a pointer to adjoint side."
//
// This test IS that arrangement: two clients, each holding ONLY the half it uses, both faces coming from
// ONE concrete integrator.  The adjointness identity still holds across them -- which is the point:
// splitting the INTERFACE does not split the OBJECT, so the pairing survives while each client's
// dependency shrinks to what it actually calls.
namespace {

//! Stands in for the CHARGE DENSITY: owns D as private state, and contracts it against a forward handle a
//! caller supplies.  It cannot reach an adjoint, and does not want one.
class DensityLike
{
public:
    explicit DensityLike(hmat_t<double> D) : itsD(std::move(D)) {}
    rvec_t ProjectOnto(const qcMesh::MatrixForward<double>& fwd) const {return fwd.Forward(itsD);}
private:
    hmat_t<double> itsD;                       // private: only the density contracts its own D
};

//! Stands in for the BASIS/term side: assembles a matrix from a field, per block.  It cannot collocate a
//! density, and does not want to.
class AssemblerLike
{
public:
    hmat_t<double> Assemble(const qcMesh::MatrixAdjoint<double>& adj, const rvec_t& v) const
        {return adj.Adjoint(v);}
    double Energy(const qcMesh::MatrixAdjoint<double>& adj, const rvec_t& f) const
        {return adj.Integrate(f);}
};

} // anon

TEST(MatrixIntegrator, TheTwoHalvesServeSeparateClientsAndStillPair)
{
    const qcMesh::Mesh m=SpreadMesh(20);
    const TinyBasis    b;
    const qcMesh::DenseMatrixIntegrator<double> I(m, b);      // ONE concrete object...

    hmat_t<double> D(3);
    D(0,0)=0.8; D(0,1)=0.3; D(0,2)=-0.5;
            D(1,1)=1.4; D(1,2)=0.2;
                    D(2,2)=0.6;

    const DensityLike   density(D);                           // ...two clients, each holding one face
    const AssemblerLike assembler;
    const qcMesh::MatrixForward<double>& fwd = I;
    const qcMesh::MatrixAdjoint<double>& adj = I;

    // The two faces must describe the SAME point set -- virtual inheritance makes NumPoints one function,
    // and this is what a client sizing an array relies on.
    EXPECT_EQ(fwd.NumPoints(), adj.NumPoints());

    rvec_t v(m.size());
    for (size_t g=0; g<m.size(); g++) v[g]=std::cos(0.9*double(g))-0.3;

    const rvec_t rho=density.ProjectOnto(fwd);
    rvec_t vrho(m.size());
    for (size_t g=0; g<m.size(); g++) vrho[g]=v[g]*rho[g];
    const double lhs=assembler.Energy(adj, vrho);

    const hmat_t<double> H=assembler.Assemble(adj, v);
    double rhs=0.0;
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) rhs+=H(i,j)*D(i,j);

    EXPECT_NEAR(lhs, rhs, 1e-12*std::max(1.0,std::abs(lhs)))
        << "the pairing did not survive the interface split -- two faces of one object must stay adjoint";
}

//---------------------------------------------------------------------------------------------------
// ★ THE FACTORED FORWARD IS NOT A CAPABILITY (user, 2026-09-09: "any implementation of MatrixForward can
// support both unfactored and factored forward.  Which one (or if both) gets used should be a non-issue").
// The base supplies a default that forms D = L L^dagger and delegates, so EVERY implementation answers
// both and the two agree exactly.  Overriding it is an optimisation, never a capability -- which is why
// there is no CanForwardFactored() to ask.
TEST(MatrixIntegrator, TheFactoredForwardAgreesWithTheUnfactoredOne)
{
    const qcMesh::Mesh m=SpreadMesh(18);
    const TinyBasis    b;
    const qcMesh::DenseMatrixIntegrator<double> I(m, b);   // does NOT override the factored overload

    // A rank-2 factor over a 3-function basis: D = L L^T is PSD by construction, like a real density.
    mat_t<double> L(3,2);
    L(0,0)= 0.7; L(0,1)=-0.2;
    L(1,0)=-1.1; L(1,1)= 0.4;
    L(2,0)= 0.3; L(2,1)= 0.9;

    hmat_t<double> D(3);
    for (size_t i=0;i<3;i++) for (size_t j=i;j<3;j++)
    { double s=0; for (size_t k=0;k<2;k++) s+=L(i,k)*L(j,k); D(i,j)=s; }

    const rvec_t rhoD=I.Forward(D);
    const rvec_t rhoL=I.Forward(L);
    ASSERT_EQ(rhoD.size(), rhoL.size());
    for (size_t g=0; g<rhoD.size(); g++)
        EXPECT_NEAR(rhoL[g], rhoD[g], 1e-13*std::max(1.0,std::abs(rhoD[g])))
            << "the factored forward must be the SAME map, not a second one";
}

