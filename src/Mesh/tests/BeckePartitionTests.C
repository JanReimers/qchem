// File: BeckePartitionTests.C  The molecular Becke partition weight and its geometry gradient (qchem.Mesh.BeckePartition).  The gradient is
// Enzyme's (AD of the weight, in place); it is checked against a central finite difference of the SAME weight with the point riding with
// its owner, against the HAND DERIVATION below (the ORACLE: test code only, never a production path -- the user's ruling 2026-11-08,
// derivatives are Enzyme's), and for translation invariance.  Without Enzyme the gradient refuses with NoGradient.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
import qchem.Mesh.BeckePartition;

using namespace qchem;

namespace
{
double& C(rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
double  C(const rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}

// THE ORACLE: the hand derivation the AD gradient replaced in qchem.Mesh.BeckePartition (lane B, 2026-11-07), kept as test code.
//   ∂|r − R_i|/∂R_D = u_i (δ_AD − δ_iD),  u_i = (r − R_i)/|r − R_i|;   ∂|R_i − R_j|/∂R_D = e_ij (δ_iD − δ_jD),  e_ij = (R_i − R_j)/|R_i − R_j|
//   ∂μ_ij/∂R_D = (∂|r − R_i| − ∂|r − R_j|)/R_ij − μ_ij ∂R_ij/R_ij, then the product and quotient rules through P and w.
//   ds/dμ = −½ Π_{i=0}^{k} f′(μ_i), μ_0 = μ, μ_{i+1} = f(μ_i).
double HandCutoffDerivative(double mu, int k)
{
    double d=1.0;
    for (int i=k; i>=0; i--) { d*=1.5*(1.0-mu*mu); mu=0.5*(3*mu - mu*mu*mu); }
    return -0.5*d;
}
std::vector<rvec3_t> HandGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k)
{
    const size_t n=R.size();
    std::vector<rvec3_t> out(n, rvec3_t(0,0,0));
    if (n==1) return out;
    std::vector<double> d(n);
    std::vector<rvec3_t> u(n);
    for (size_t i=0;i<n;i++) { d[i]=norm(r-R[i]); u[i] = d[i]>0.0 ? (1.0/d[i])*(r-R[i]) : rvec3_t(0,0,0); }
    std::vector<double> s(n*n, 1.0), sp(n*n, 0.0);
    std::vector<rvec3_t> dmu(n*n*n, rvec3_t(0,0,0));   // dmu[(i*n+j)*n+D]
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) if (j!=i)
    {
        const double Rij=norm(R[i]-R[j]);
        if (!(Rij>0.0)) { s[i*n+j]=qcMesh::BeckeCutoff(0.0,k); continue; }
        const rvec3_t e=(1.0/Rij)*(R[i]-R[j]);
        const double mu=(d[i]-d[j])/Rij;
        s[i*n+j]=qcMesh::BeckeCutoff(mu,k); sp[i*n+j]=HandCutoffDerivative(mu,k);
        for (size_t D=0;D<n;D++)
        {
            const double a=(owner==D)-(i==D), b=(owner==D)-(j==D), c=double(i==D)-double(j==D);
            dmu[(i*n+j)*n+D]=(1.0/Rij)*(a*u[i]-b*u[j]) - (mu*c/Rij)*e;
        }
    }
    std::vector<double> P(n, 1.0);
    std::vector<rvec3_t> dP(n*n, rvec3_t(0,0,0));
    for (size_t i=0;i<n;i++)
    {
        for (size_t j=0;j<n;j++) if (j!=i) P[i]*=s[i*n+j];
        for (size_t j=0;j<n;j++) if (j!=i)
        {
            double others=1.0;
            for (size_t m=0;m<n;m++) if (m!=i && m!=j) others*=s[i*n+m];
            const double f=sp[i*n+j]*others;
            if (f!=0.0) for (size_t D=0;D<n;D++) dP[i*n+D]+=f*dmu[(i*n+j)*n+D];
        }
    }
    double Z=0.0; for (double p : P) Z+=p;
    if (!(Z>0.0)) return out;
    for (size_t D=0;D<n;D++)
    {
        rvec3_t dZ(0,0,0); for (size_t i=0;i<n;i++) dZ+=dP[i*n+D];
        out[D]=(1.0/Z)*dP[owner*n+D] - (P[owner]/(Z*Z))*dZ;
    }
    return out;
}
// A bent triatomic plus a fourth atom off its plane: no symmetry to hide a sign error behind
const std::vector<rvec3_t> kAtoms={rvec3_t(0.0,0.0,0.0), rvec3_t(1.43,0.0,1.11), rvec3_t(-1.43,0.0,1.11), rvec3_t(0.4,2.1,-0.7)};
const std::vector<rvec3_t> kOffsets={rvec3_t(0.3,-0.2,0.1), rvec3_t(-0.9,0.6,0.8), rvec3_t(1.4,1.1,-0.5), rvec3_t(0.05,0.02,-0.04), rvec3_t(-2.0,0.3,1.7)};
}

TEST(BeckePartition, TheOraclesCutoffDerivativeIsTheIteratedChainRule)
{
    for (int k=0;k<=3;++k)
        for (double mu : {-0.95,-0.6,-0.2,0.0,0.15,0.5,0.9})
        {
            const double h=1e-6, fd=(qcMesh::BeckeCutoff(mu+h,k)-qcMesh::BeckeCutoff(mu-h,k))/(2*h);
            EXPECT_NEAR(HandCutoffDerivative(mu,k), fd, 1e-8) << "k " << k << " mu " << mu;
        }
}

TEST(BeckePartition, TheWeightsAreAPartitionOfUnity)
{
    for (const auto& off : kOffsets)
        for (int k : {2,3})
        {
            const rvec3_t r=kAtoms[1]+off;
            double sum=0.0;
            for (size_t A=0;A<kAtoms.size();++A) sum+=qcMesh::BeckePartitionWeight(kAtoms, A, r, k);
            EXPECT_NEAR(sum, 1.0, 1e-14);
        }
    EXPECT_DOUBLE_EQ(qcMesh::BeckePartitionWeight({rvec3_t(0,0,0)}, 0, rvec3_t(1,2,3), 2), 1.0);
}

TEST(BeckePartition, TheGradientIsTheFiniteDifferenceAndTheOracleWithThePointRidingWithItsOwner)
{
    const double h=1e-5;
    double largest=0.0;   // the points must sit where the partition actually moves, or the test proves nothing
    for (int k : {2,3})
        for (size_t A=0;A<kAtoms.size();++A)
            for (const auto& off : kOffsets)
            {
                const auto g=qcMesh::BeckePartitionGradient(kAtoms, A, kAtoms[A]+off, k);
                const auto hand=HandGradient(kAtoms, A, kAtoms[A]+off, k);
                ASSERT_EQ(g.size(), kAtoms.size());
                rvec3_t sum(0,0,0);
                for (size_t D=0;D<kAtoms.size();++D)
                {
                    sum+=g[D]; largest=std::max(largest, norm(g[D]));
                    EXPECT_LT(norm(g[D]-hand[D]), 1e-13) << "the AD gradient is the oracle's: k " << k << " owner " << A << " atom " << D;
                    for (size_t x=0;x<3;++x)
                    {
                        auto Rp=kAtoms, Rm=kAtoms; C(Rp[D],x)+=h; C(Rm[D],x)-=h;
                        const double fd=(qcMesh::BeckePartitionWeight(Rp, A, Rp[A]+off, k)-qcMesh::BeckePartitionWeight(Rm, A, Rm[A]+off, k))/(2*h);
                        EXPECT_NEAR(C(g[D],x), fd, 1e-8) << "k " << k << " owner " << A << " atom " << D << " axis " << x;
                    }
                }
                EXPECT_LT(norm(sum), 1e-13) << "translation invariance: the gradient sums to zero over the atoms";
            }
    EXPECT_GT(largest, 0.1) << "the sampled points are in the partition's transition regions";
    const auto one=qcMesh::BeckePartitionGradient({rvec3_t(0,0,0)}, 0, rvec3_t(1,0,0), 2);
    ASSERT_EQ(one.size(), 1u);
    EXPECT_EQ(norm(one[0]), 0.0) << "one atom: the weight is 1 everywhere";
    const auto atNucleus=qcMesh::BeckePartitionGradient(kAtoms, 1, kAtoms[1], 3);
    for (size_t D=0;D<kAtoms.size();++D) EXPECT_LT(norm(atNucleus[D]-HandGradient(kAtoms, 1, kAtoms[1], 3)[D]), 1e-15) << "at the owner's nucleus: zero";
}

// THE PERIODIC PARTITION's GRADIENT (lane B 2026-11-19): PeriodicBeckeWeightGradient against the central difference of PeriodicBeckeWeight
// itself, the point riding its owner (r = R_ia + u), on an oblique cell with two atoms -- per CELL ATOM (every image of an atom moves with
// it: the force), and per homogeneous strain (cell, atoms and images strained, the template offset u held: the virial Σ_X g_X ⊗ X).  The
// eps screen is set at 1e-12 so no P-set member crosses it inside ±h.
namespace
{
//! A general cell: columns of A are the lattice vectors (ToCartesian(f) = A f).
struct TestCell
{
    double A[3][3];
    rvec3_t ToCartesian(const rvec3_t& f) const
    {
        return rvec3_t(A[0][0]*f.x+A[0][1]*f.y+A[0][2]*f.z, A[1][0]*f.x+A[1][1]*f.y+A[1][2]*f.z, A[2][0]*f.x+A[2][1]*f.y+A[2][2]*f.z);
    }
    rvec3_t ToFractional(const rvec3_t& r) const
    {
        const double (&a)[3][3]=A;
        const double det=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
        double inv[3][3];
        inv[0][0]=(a[1][1]*a[2][2]-a[1][2]*a[2][1])/det; inv[0][1]=(a[0][2]*a[2][1]-a[0][1]*a[2][2])/det; inv[0][2]=(a[0][1]*a[1][2]-a[0][2]*a[1][1])/det;
        inv[1][0]=(a[1][2]*a[2][0]-a[1][0]*a[2][2])/det; inv[1][1]=(a[0][0]*a[2][2]-a[0][2]*a[2][0])/det; inv[1][2]=(a[0][2]*a[1][0]-a[0][0]*a[1][2])/det;
        inv[2][0]=(a[1][0]*a[2][1]-a[1][1]*a[2][0])/det; inv[2][1]=(a[0][1]*a[2][0]-a[0][0]*a[2][1])/det; inv[2][2]=(a[0][0]*a[1][1]-a[0][1]*a[1][0])/det;
        return rvec3_t(inv[0][0]*r.x+inv[0][1]*r.y+inv[0][2]*r.z, inv[1][0]*r.x+inv[1][1]*r.y+inv[1][2]*r.z, inv[2][0]*r.x+inv[2][1]*r.y+inv[2][2]*r.z);
    }
    //! The minimum interplanar spacing 1/max_i |b_i|, b_i the reciprocal rows (no 2π).
    double DPlane() const
    {
        double m=0.0;   // the rows of A^-1, probed column by column: row i = ((A^-1 e_0)_i, (A^-1 e_1)_i, (A^-1 e_2)_i)
        const rvec3_t c0=ToFractional(rvec3_t(1,0,0)), c1=ToFractional(rvec3_t(0,1,0)), c2=ToFractional(rvec3_t(0,0,1));
        for (int i=0;i<3;i++)
        {
            const rvec3_t row(C(c0,i), C(c1,i), C(c2,i));
            m=std::max(m, norm(row));
        }
        return 1.0/m;
    }
};
TestCell Strained(const TestCell& c, int u, int v, double s)   // (1 + s e_uv) A
{
    TestCell o=c;
    for (int j=0;j<3;j++) o.A[u][j]=c.A[u][j]+s*c.A[v][j];
    return o;
}
} // namespace

TEST(BeckePartition, ThePeriodicWeightsGradientIsItsFiniteDifferencePerAtomAndPerStrain)
{
    const TestCell cell{{{6.0,0.9,0.4},{0.3,5.5,0.7},{-0.2,0.5,6.4}}};
    const rvec3_t f[2]={rvec3_t(0.12,0.21,0.08), rvec3_t(0.47,0.55,0.61)};
    auto sites=[&](const TestCell& c){ std::vector<rvec3_t> R; for (const auto& fi : f) R.push_back(c.ToCartesian(fi)); return R; };
    const int k=3; const double eps=1e-12;
    const size_t ia=1;
    const rvec3_t u(0.9,-0.6,1.3);   // a template offset, far enough out that the partition is not flat
    auto weight=[&](const TestCell& c, const std::vector<rvec3_t>& R)
    {
        std::vector<qcMesh::BeckeImage> im; long t=-1; qcMesh::BeckePointCost pc;
        return qcMesh::PeriodicBeckeWeight(c, R, c.DPlane(), ia, R[ia]+u, k, eps, im, t, pc, false);
    };
    const std::vector<rvec3_t> R0=sites(cell);
    std::vector<qcMesh::BeckeImage> im; long t=-1; qcMesh::BeckePointCost pc;
    const double w0=qcMesh::PeriodicBeckeWeight(cell, R0, cell.DPlane(), ia, R0[ia]+u, k, eps, im, t, pc, false);
    ASSERT_GT(w0, 1e-3); ASSERT_LT(w0, 1.0-1e-3);
    ASSERT_GE(t, 0);
    const std::vector<rvec3_t> g=qcMesh::PeriodicBeckeWeightGradient(im, size_t(t), u, k);
    ASSERT_EQ(g.size(), im.size());
    std::vector<rvec3_t> perAtom(2, rvec3_t(0,0,0));
    double V[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    rvec3_t sum(0,0,0);
    for (size_t i=0;i<im.size();i++)
    {
        perAtom[im[i].atom]+=g[i]; sum+=g[i];
        for (size_t a=0;a<3;a++) for (size_t b=0;b<3;b++) V[a][b]+=C(g[i],a)*C(im[i].R,b);
    }
    EXPECT_LT(norm(sum), 1e-13) << "translation: the image gradients sum to zero";
    const double h=1e-5;
    double big=0.0, worst=0.0;
    for (size_t A=0;A<2;A++) for (size_t x=0;x<3;x++)
    {
        auto at=[&](double s){ std::vector<rvec3_t> R=R0; C(R[A],x)+=s; return weight(cell, R); };
        const double fd=(at(h)-at(-h))/(2*h), an=C(perAtom[A],x);
        big=std::max(big,std::fabs(fd)); worst=std::max(worst,std::fabs(an-fd));
        EXPECT_NEAR(an, fd, 1e-8) << "atom " << A << " axis " << x;
    }
    for (int a=0;a<3;a++) for (int b=0;b<3;b++)
    {
        auto at=[&](double s){ const TestCell c=Strained(cell,a,b,s); return weight(c, sites(c)); };
        const double fd=(at(h)-at(-h))/(2*h);
        big=std::max(big,std::fabs(fd)); worst=std::max(worst,std::fabs(V[a][b]-fd));
        EXPECT_NEAR(V[a][b], fd, 1e-8) << "strain " << a << b;
    }
    std::printf("  [periodic Becke] w %.6f  %zu images  max|dw| %.3e  worst |analytic - FD| %.2e\n", w0, im.size(), big, worst);
    EXPECT_GT(big, 1e-3) << "a flat partition cannot show the gradient is wired";
}
