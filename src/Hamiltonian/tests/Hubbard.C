// File: Hubbard.C  DFT+U unit gates -- the LowdinProjector pair on a synthetic block, no SCF.
//
// The property the term rests on is that the occupation-matrix FORWARD and the potential ADJOINT are the
// same T = S^{1/2}[:,M] contracted both ways, so H = dE/dD holds exactly; and that Loewdin occupations of
// an idempotent density (projected onto a manifold that spans it) are exactly 0 or 1.
//
// INCREMENT 2 (the ManifoldSymmetry gates): a d shell under O_h splits {2, 3} = e_g + t2g with the textbook
// characters; under the D_3d site group of an AFM-II MnO Mn inside grey O_h it splits {1, 2, 2} =
// a1g < t2g, e_g < e_g, e_g < t2g; a (U_eg, U_t2g) Uirrep reproduces E_U and W by hand; equal Uirrep is
// bit-identical with the shell-averaged U; and the seed n = 0 (one big degenerate cluster) is still named
// with purity 1 -- the LAPACK basis is rotated inside the cluster, never trusted.
#include "gtest/gtest.h"
#include <algorithm>
#include <complex>
#include <cmath>
#include <memory>
#include <vector>
import qchem.Hamiltonian.Internal.Hubbard;   // LowdinProjector, ManifoldSymmetry (unit tests may import Internal)
import qchem.Symmetry.Molecule.SphericalRep;  // SphericalShellRep (a d shell's operation rep)
import qchem.Symmetry.Molecule.OperationRep;  // AoShell
import qchem.Math.Angular;                    // Math::SphericalShell(2)
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::Hamiltonian;

namespace
{
// A 4-function block whose overlap is a mildly non-orthogonal SPD matrix; manifold = functions {1,2}.
hmat_t<double> OverlapR()
{
    hmat_t<double> S(4);
    for (size_t i=0;i<4;i++) S(i,i)=1.0;
    S(0,1)=0.2; S(0,2)=-0.1; S(0,3)=0.05; S(1,2)=0.3; S(1,3)=-0.15; S(2,3)=0.1;
    return S;
}
hmat_t<dcmplx> OverlapC()
{
    hmat_t<dcmplx> S(4);
    for (size_t i=0;i<4;i++) S(i,i)=1.0;
    S(0,1)=dcmplx(0.2,0.05); S(0,2)=dcmplx(-0.1,0.0); S(0,3)=dcmplx(0.05,-0.02);
    S(1,2)=dcmplx(0.3,0.1);  S(1,3)=dcmplx(-0.15,0.0); S(2,3)=dcmplx(0.1,0.03);
    return S;
}
const std::vector<std::vector<size_t>> M{{1,2}};
}

//---------------------------------------------------------------------------------------------------
TEST(LowdinProjector, ForwardAndAdjointAreExactlyAdjoint_Real)
{
    const LowdinProjector<double> P(OverlapR(), M);
    EXPECT_EQ(P.NumCoefficients(), 4u);   // 2x2 flattened
    hmat_t<double> D(4); D(0,0)=1.3; D(0,1)=-0.7; D(0,2)=0.2; D(0,3)=0.1; D(1,1)=2.1; D(1,2)=0.4; D(1,3)=-0.3; D(2,2)=0.6; D(2,3)=0.05; D(3,3)=0.9;
    const rvec_t W{0.9, -0.4, -0.4, 1.7};   // a symmetric 2x2, flattened
    const rvec_t n=P.Forward(D);
    double lhs=0.0; for (size_t a=0;a<4;a++) lhs+=W[a]*n[a];
    const hmat_t<double> V=P.Adjoint(W);
    double rhs=0.0; for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++) rhs+=V(i,j)*D(i,j);
    EXPECT_NEAR(lhs, rhs, 1e-13) << "Tr(W n) != Tr(V D): the pair is not adjoint";
}

TEST(LowdinProjector, ForwardAndAdjointAreExactlyAdjoint_Complex)
{
    const LowdinProjector<dcmplx> P(OverlapC(), M);
    hmat_t<dcmplx> D(4);
    D(0,0)=1.3; D(0,1)=dcmplx(-0.7,0.2); D(0,2)=dcmplx(0.2,-0.1); D(0,3)=0.1; D(1,1)=2.1; D(1,2)=dcmplx(0.4,0.3); D(1,3)=-0.3; D(2,2)=0.6; D(2,3)=dcmplx(0.05,0.02); D(3,3)=0.9;
    const rvec_t W{0.9, -0.4, -0.4, 1.7};
    const rvec_t n=P.Forward(D);
    double lhs=0.0; for (size_t a=0;a<4;a++) lhs+=W[a]*n[a];
    const hmat_t<dcmplx> V=P.Adjoint(W);
    dcmplx rhs=0.0; for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++) rhs+=V(i,j)*std::conj(dcmplx(D(i,j)));
    // With a REAL symmetric W the forward keeps only Re(T^dag D T), so the adjoint identity holds on the
    // real part; Tr(V D) is real for Hermitian V, D.
    EXPECT_NEAR(lhs, std::real(rhs), 1e-13);
    EXPECT_NEAR(std::imag(rhs), 0.0, 1e-13);
}

//---------------------------------------------------------------------------------------------------
// An idempotent density living entirely inside the manifold: D = S^{-1/2} P S^{-1/2} with P a projector
// onto Loewdin-orthogonal functions {1}. Its Loewdin occupation matrix is P restricted -- eigenvalues {1,0}.
TEST(LowdinProjector, AnIdempotentDensityInTheManifoldHasOccupations0And1)
{
    const hmat_t<double> S=OverlapR();
    rvec_t w; rmat_t V; blazem::eigen(S, w, V);
    rmat_t Sinvh(4,4);
    for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++)
    { double s=0; for (size_t k=0;k<4;k++) s+=V(i,k)*V(j,k)/std::sqrt(w[k]); Sinvh(i,j)=s; }
    // D = S^{-1/2} e1 e1^T S^{-1/2}: one electron in Loewdin function 1 (a member of the manifold).
    hmat_t<double> D(4);
    for (size_t i=0;i<4;i++) for (size_t j=i;j<4;j++) D(i,j)=Sinvh(i,1)*Sinvh(j,1);
    const LowdinProjector<double> P(S, M);
    const rvec_t n=P.Forward(D);           // the 2x2 block over {1,2}: should be [[1,0],[0,0]]
    EXPECT_NEAR(n[0], 1.0, 1e-12);
    EXPECT_NEAR(n[1], 0.0, 1e-12);
    EXPECT_NEAR(n[2], 0.0, 1e-12);
    EXPECT_NEAR(n[3], 0.0, 1e-12);
    EXPECT_NEAR(P.Integrate(n), 1.0, 1e-12);   // the manifold charge
}

//---------------------------------------------------------------------------------------------------
TEST(LowdinProjector, TheFactoredForwardAgreesWithTheUnfactoredOne)
{
    const LowdinProjector<double> P(OverlapR(), M);
    rmat_t L(4,2); L(0,0)=1.0; L(0,1)=0.2; L(1,0)=-0.5; L(1,1)=1.5; L(2,0)=0.3; L(2,1)=-0.4; L(3,0)=0.1; L(3,1)=0.7;
    hmat_t<double> D(4);
    for (size_t i=0;i<4;i++) for (size_t j=i;j<4;j++) { double s=0; for (size_t m=0;m<2;m++) s+=L(i,m)*L(j,m); D(i,j)=s; }
    const rvec_t a=P.Forward(D), b=P.Forward(L);
    ASSERT_EQ(a.size(), b.size());
    for (size_t k=0;k<a.size();k++) EXPECT_NEAR(a[k], b[k], 1e-13);
}

//=============================================================================================== increment 2
namespace
{
using Symmetry::Molecule::AoShell;
using Symmetry::Molecule::SphericalShellRep;

//! ONE d shell as a basis reports it: the raw real harmonics' rep, and per-component norms N = 1/|chi_m|
//! (angular norms^2 of xy, yz, z^2-(x^2+y^2)/2, xz, x^2-y^2 are 4pi/15, 4pi/15, 4pi/5, 4pi/15, 16pi/15 --
//! the ratio 1 : 3 : 4 that makes the RAW rep non-orthogonal and the normalised one orthogonal).
AoShell DShell(size_t offset=0)
{
    AoShell sh;
    sh.shellType=2; sh.center=rvec3_t(0,0,0); sh.offset=offset;
    sh.rep=std::make_shared<SphericalShellRep>(Math::SphericalShell(2));
    sh.norm=rvec_t{1.0, 1.0, 1.0/std::sqrt(3.0), 1.0, 0.5};
    return sh;
}

//! O_h as the 48 signed permutation matrices.
std::vector<rmat3d_t> Oh()
{
    std::vector<rmat3d_t> ops;
    int perm[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (auto& pm : perm)
        for (int sgn=0; sgn<8; sgn++)
        {
            rmat3d_t R(0,0,0, 0,0,0, 0,0,0);
            for (int r=0;r<3;r++) R(r+1, pm[r]+1) = ((sgn>>r)&1) ? -1.0 : 1.0;   // Matrix3D is 1-based
            ops.push_back(R);
        }
    return ops;
}
//! D_3d as the subgroup of O_h fixing the [111] axis (up to sign): 12 ops.
std::vector<rmat3d_t> D3d()
{
    std::vector<rmat3d_t> ops;
    for (const rmat3d_t& R : Oh())
    {
        const rvec3_t a(1,1,1), Ra=R*a;
        if (norm(Ra-a)<1e-12 || norm(Ra+a)<1e-12) ops.push_back(R);
    }
    return ops;
}
//! A symmetric n = Sum_k c_k P_k over the site projectors (exactly site-symmetric, degenerate per irrep).
rsmat_t FromProjectors(const ManifoldSymmetry& sym, const std::vector<double>& c)
{
    const size_t m=sym.Size();
    rmat_t n(m,m,0.0);
    for (size_t k=0;k<c.size();k++) n+=c[k]*sym.SiteProjector(k);
    rsmat_t out(m);
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++) out(a,b)=n(a,b);
    return out;
}
double Trace(const rmat_t& A) { double t=0; for (size_t i=0;i<A.rows();i++) t+=A(i,i); return t; }
double MaxAbs(const rmat_t& A) { double m=0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.columns();j++) m=std::max(m, std::abs(A(i,j))); return m; }
rmat_t Identity(size_t m) { rmat_t I(m,m,0.0); for (size_t i=0;i<m;i++) I(i,i)=1.0; return I; }
}

TEST(ManifoldSymmetry, ADShellUnderOhSplitsEgPlusT2gWithTheTextbookCharacters)
{
    const std::vector<rmat3d_t> ops=Oh();
    ASSERT_EQ(ops.size(), 48u);
    const std::vector<rmat_t> D=ManifoldSymmetry::Rep({DShell()}, ops);
    for (const rmat_t& d : D) EXPECT_NEAR(MaxAbs(blazem::trans(d)*d-Identity(5)), 0.0, 1e-12);   // orthogonal in the normalised basis
    const ManifoldSymmetry sym(D, D);
    ASSERT_EQ(sym.Size(), 5u);
    ASSERT_EQ(sym.Irreps().size(), 2u);                     // e_g, t2g -- sorted by dimension
    EXPECT_EQ(sym.Irreps()[0].dim, 2u);
    EXPECT_EQ(sym.Irreps()[1].dim, 3u);
    // Characters on three named classes: C_4 about z (t2g -1, e_g 0), C_3 about [111] (t2g 0, e_g -1),
    // inversion (both g: +3, +2), C_2 about z (t2g -1, e_g +2).
    auto same=[](const rmat3d_t& A, const rmat3d_t& B){ for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::abs(A(i,j)-B(i,j))>1e-12) return false; return true; };
    auto find=[&](const rmat3d_t& R){ for (size_t g=0;g<ops.size();g++) if (same(ops[g],R)) return g; throw std::logic_error("op"); };
    const size_t c4=find(rmat3d_t(0,-1,0, 1,0,0, 0,0,1)), c3=find(rmat3d_t(0,0,1, 1,0,0, 0,1,0)),
                 inv=find(rmat3d_t(-1,0,0, 0,-1,0, 0,0,-1)), c2=find(rmat3d_t(-1,0,0, 0,-1,0, 0,0,1));
    const rvec_t& eg=sym.Irreps()[0].chi; const rvec_t& t2g=sym.Irreps()[1].chi;
    EXPECT_NEAR(eg[c4],  0.0, 1e-9);  EXPECT_NEAR(t2g[c4], -1.0, 1e-9);
    EXPECT_NEAR(eg[c3], -1.0, 1e-9);  EXPECT_NEAR(t2g[c3],  0.0, 1e-9);
    EXPECT_NEAR(eg[inv], 2.0, 1e-9);  EXPECT_NEAR(t2g[inv], 3.0, 1e-9);
    EXPECT_NEAR(eg[c2],  2.0, 1e-9);  EXPECT_NEAR(t2g[c2], -1.0, 1e-9);
    // Same group twice: the slots are the irreps themselves, {2, 3}, ungraded.
    EXPECT_FALSE(sym.Graded());
    ASSERT_EQ(sym.Slots().size(), 2u);
    EXPECT_EQ(sym.Slots()[0].dim, 2u); EXPECT_EQ(sym.Slots()[0].irrep, 0u);
    EXPECT_EQ(sym.Slots()[1].dim, 3u); EXPECT_EQ(sym.Slots()[1].irrep, 1u);
    // The projectors are projectors, orthogonal, and sum to the identity.
    for (size_t k=0;k<2;k++)
    {
        const rmat_t& P=sym.SiteProjector(k);
        EXPECT_NEAR(MaxAbs(P*P-P), 0.0, 1e-10);
        EXPECT_NEAR(Trace(P), double(sym.Irreps()[k].dim), 1e-10);
    }
    EXPECT_NEAR(MaxAbs(sym.SiteProjector(0)*sym.SiteProjector(1)), 0.0, 1e-10);
    EXPECT_NEAR(MaxAbs(sym.SiteProjector(0)+sym.SiteProjector(1)-Identity(5)), 0.0, 1e-10);
}

TEST(ManifoldSymmetry, AnOhSymmetricOccupationGetsATwoSlotUirrepEnergyAndPotentialByHand)
{
    const std::vector<rmat_t> D=ManifoldSymmetry::Rep({DShell()}, Oh());
    const ManifoldSymmetry sym(D, D);
    const double b=0.35, a=0.90;                             // e_g and t2g occupations (slot order: e_g first)
    const rsmat_t n=FromProjectors(sym, {b, a});
    rvec_t lam; rmat_t v;
    const auto L=sym.Label(n, lam, v);
    ASSERT_EQ(L.slot.size(), 5u);
    EXPECT_NEAR(L.purity, 1.0, 1e-10);
    rvec_t byslot(2, 0.0); std::vector<int> count(2,0);
    for (size_t i=0;i<5;i++) { byslot[L.slot[i]]+=lam[i]; count[L.slot[i]]++; }
    EXPECT_EQ(count[0], 2); EXPECT_EQ(count[1], 3);
    EXPECT_NEAR(byslot[0], 2*b, 1e-12); EXPECT_NEAR(byslot[1], 3*a, 1e-12);
    // Hand values: E_U = 2 U_eg/2 b(1-b) + 3 U_t2g/2 a(1-a);  W = U_eg(1/2-b) P_eg + U_t2g(1/2-a) P_t2g.
    const double Ueg=0.10, Ut2g=0.20;
    rmat_t W;
    const double EU=DudarevInEigenbasis(lam, v, L.slot, {Ueg, Ut2g}, 0.0, W);
    EXPECT_NEAR(EU, Ueg*b*(1-b) + 1.5*Ut2g*a*(1-a), 1e-13);
    const rmat_t Whand=Ueg*(0.5-b)*sym.SiteProjector(0)+Ut2g*(0.5-a)*sym.SiteProjector(1);
    EXPECT_NEAR(MaxAbs(W-Whand), 0.0, 1e-12);
    // Equal Uirrep IS the shell-averaged U, bit for bit.
    rmat_t W1, W2;
    const double E1=DudarevInEigenbasis(lam, v, L.slot, {0.147, 0.147}, 0.0, W1);
    const double E2=DudarevInEigenbasis(lam, v, {},     {},             0.147, W2);
    EXPECT_EQ(E1, E2);
    EXPECT_EQ(MaxAbs(W1-W2), 0.0);
}

TEST(ManifoldSymmetry, ADShellOnAD3dSiteInsideGreyOhHasThreeSlots_1_2_2_WithTheirParents)
{
    const std::vector<rmat3d_t> site=D3d();
    ASSERT_EQ(site.size(), 12u);
    const AoShell d=DShell();
    const ManifoldSymmetry sym(ManifoldSymmetry::Rep({d}, site), ManifoldSymmetry::Rep({d}, Oh()));
    EXPECT_TRUE(sym.Graded());
    ASSERT_EQ(sym.Irreps().size(), 2u);                     // a1g (1), e_g (2, twice)
    EXPECT_EQ(sym.Irreps()[0].dim, 1u);
    EXPECT_EQ(sym.Irreps()[1].dim, 2u);
    ASSERT_EQ(sym.Grey().size(), 2u);                       // e_g (2), t2g (3)
    ASSERT_EQ(sym.Slots().size(), 3u);
    // sorted by (dim, parent, irrep): [0] a1g < t2g, [1] e_g < e_g, [2] e_g < t2g
    EXPECT_EQ(sym.Slots()[0].dim, 1u); EXPECT_EQ(sym.Slots()[0].irrep, 0u); EXPECT_EQ(sym.Slots()[0].parent, 1u);
    EXPECT_EQ(sym.Slots()[1].dim, 2u); EXPECT_EQ(sym.Slots()[1].irrep, 1u); EXPECT_EQ(sym.Slots()[1].parent, 0u);
    EXPECT_EQ(sym.Slots()[2].dim, 2u); EXPECT_EQ(sym.Slots()[2].irrep, 1u); EXPECT_EQ(sym.Slots()[2].parent, 1u);
    // The a1g function is d_{z^2} along [111]: its projector is rank 1 and lies inside t2g.
    const rmat_t& Pa=sym.SiteProjector(0);
    EXPECT_NEAR(Trace(Pa), 1.0, 1e-10);
    EXPECT_NEAR(MaxAbs(sym.GreyProjector(1)*Pa-Pa), 0.0, 1e-10);
    // Eight d shells (MnO's VA span): the slot dimensions scale, the table does not change shape.
    std::vector<AoShell> eight; for (size_t k=0;k<8;k++) eight.push_back(DShell(5*k));
    const ManifoldSymmetry sym8(ManifoldSymmetry::Rep(eight, site), ManifoldSymmetry::Rep(eight, Oh()));
    ASSERT_EQ(sym8.Slots().size(), 3u);
    EXPECT_EQ(sym8.Slots()[0].dim, 8u); EXPECT_EQ(sym8.Slots()[1].dim, 16u); EXPECT_EQ(sym8.Slots()[2].dim, 16u);
}

TEST(ManifoldSymmetry, TheSeedNIsOneDegenerateClusterAndIsStillNamedWithPurityOne)
{
    // n = 0: LAPACK hands back the raw harmonics, which are NOT D_3d-adapted (d_xy mixes a1g and e_g along
    // [111]).  The labelling must rotate inside the cluster and come out pure, with the {1,2,2} counts.
    const AoShell d=DShell();
    const ManifoldSymmetry sym(ManifoldSymmetry::Rep({d}, D3d()), ManifoldSymmetry::Rep({d}, Oh()));
    rsmat_t n(5); for (size_t a=0;a<5;a++) for (size_t b=a;b<5;b++) n(a,b)=0.0;
    rvec_t lam; rmat_t v;
    const auto L=sym.Label(n, lam, v);
    EXPECT_NEAR(L.purity, 1.0, 1e-10);
    std::vector<int> count(3,0); for (size_t s : L.slot) count[s]++;
    EXPECT_EQ(count[0], 1); EXPECT_EQ(count[1], 2); EXPECT_EQ(count[2], 2);
    for (double w : L.parentage) EXPECT_NEAR(w, 1.0, 1e-10);
    // The eigenvectors still diagonalise n (trivially) and stay orthonormal after the rotation.
    EXPECT_NEAR(MaxAbs(blazem::trans(v)*v-Identity(5)), 0.0, 1e-12);
    // A trigonal n that MIXES the two e_g copies (the physical case): purity stays 1 (it is D_3d-symmetric)
    // but the parentage of the e_g slots drops below 100 % -- the diagnostic the MnO trace is read by.
    const rmat_t X=sym.SiteProjector(1)*(sym.GreyProjector(0)-sym.GreyProjector(1))*sym.SiteProjector(1);   // e_g x e_g, mixes parents
    rmat_t nm=0.5*sym.SiteProjector(0)+0.6*sym.SiteProjector(1)+0.1*(X*X);   // X*X is symmetric and site-invariant
    // add an off-parent coupling: (P_eg P_t2g-parent P_eg) itself is not the identity on e_g, so nm's
    // e_g eigenvectors are parent-mixed
    nm+=0.05*(sym.SiteProjector(1)*sym.GreyProjector(1)*sym.SiteProjector(1));
    rsmat_t nS(5); for (size_t a=0;a<5;a++) for (size_t b=a;b<5;b++) nS(a,b)=0.5*(nm(a,b)+nm(b,a));
    const auto L2=sym.Label(nS, lam, v);
    EXPECT_NEAR(L2.purity, 1.0, 1e-10);
    std::vector<int> count2(3,0); for (size_t s : L2.slot) count2[s]++;
    EXPECT_EQ(count2[0], 1); EXPECT_EQ(count2[1]+count2[2], 4);
}

//=============================================================================================== increment 3
// THE CONTRACTED MANIFOLD (slice C): a manifold that is a fixed combination chi = phi[:,cols] V of the block's
// columns, projected ATOMICALLY (T = S[:,cols] Vt with Vt S-orthonormal: T^dagger c = <chi|psi>), beside the
// column manifold's Löwdin S^{1/2}.  Claims: the chi's are S-orthonormal; T^dagger c == <chi|psi> == the
// coefficient map; a density made of ONE normalised chi has occupation exactly 1 on it and 0 on an
// S-orthogonal partner; the pair is adjoint; an explicit identity contraction is Löwdin-within-the-manifold.
TEST(LowdinProjector, AContractedManifoldProjectsAtomically)
{
    // A 6-function block: two "shells" of 3 components each (columns 0-2 and 3-5), overlapping strongly.
    hmat_t<double> S(6);
    for (size_t i=0;i<6;i++) S(i,i)=1.0;
    for (size_t k=0;k<3;k++) { S(k,k+3)=0.7; }                       // shell 0 <-> shell 1, same component
    S(0,1)=0.1; S(3,4)=0.1; S(1,5)=0.05;                              // some cross terms
    const std::vector<std::vector<size_t>> cols{{0,1,2,3,4,5}};
    // The contraction: chi_m = r0 phi_{0,m} + r1 phi_{1,m}, m = 0..2.
    mat_t<double> V(6,3,0.0);
    for (size_t m=0;m<3;m++) { V(m,m)=0.6; V(3+m,m)=0.5; }
    const LowdinProjector<double> P(S, cols, {V});
    ASSERT_EQ(P.NumManifolds(), 1u);
    ASSERT_EQ(P.Size(0), 3u);
    EXPECT_EQ(P.NumCoefficients(), 9u);
    // Vt is S-orthonormal on the columns.
    const rmat_t& Vt=P.Contraction(0);
    rmat_t Sc(6,6); for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++) Sc(i,j)=S(i,j);
    const rmat_t G=blazem::trans(Vt)*Sc*Vt;
    for (size_t a=0;a<3;a++) for (size_t b=0;b<3;b++) EXPECT_NEAR(G(a,b), a==b?1.0:0.0, 1e-12);
    // T^dagger c is <chi|psi> = Vt^T S c, and Coefficients() is the same map.
    vec_t<double> c(6); for (size_t i=0;i<6;i++) c[i]=0.1*(i+1);
    const vec_t<double> l1=blazem::trans(P.T(0))*c, l2=blazem::trans(Vt)*(Sc*c), l3=P.Coefficients(0)*c;
    for (size_t a=0;a<3;a++) { EXPECT_NEAR(l1[a], l2[a], 1e-13); EXPECT_NEAR(l3[a], l2[a], 1e-13); }
    // A density that IS one normalised chi: occupation 1 on that chi, 0 on the others.
    const vec_t<double> chi0=blazem::column(Vt,0);                   // AO coefficients of chi_0 (already S-normalised)
    hmat_t<double> D(6); for (size_t i=0;i<6;i++) for (size_t j=i;j<6;j++) D(i,j)=chi0[i]*chi0[j];
    const rvec_t n=P.Forward(D);
    EXPECT_NEAR(n[0*3+0], 1.0, 1e-12); EXPECT_NEAR(n[1*3+1], 0.0, 1e-12); EXPECT_NEAR(n[2*3+2], 0.0, 1e-12);
    EXPECT_NEAR(P.Integrate(n), 1.0, 1e-12);
    // An EXPLICIT identity contraction is NOT the column manifold: it is Löwdin WITHIN the manifold (chi = phi
    // S_cc^{-1/2}, projected atomically), where the column manifold is Löwdin over the WHOLE block (S^{1/2}[:,M]).
    // Both are orthonormal sets, so the charge of a density inside the span is the same -- the occupation
    // MATRICES differ (different functions), and a caller who wants CP2K's convention passes NO contraction.
    mat_t<double> I6(6,6,0.0); for (size_t i=0;i<6;i++) I6(i,i)=1.0;
    const LowdinProjector<double> Pc(S, cols), Pi(S, cols, {I6});
    const rvec_t nc=Pc.Forward(D), ni=Pi.Forward(D);
    EXPECT_NEAR(Pc.Integrate(nc), 1.0, 1e-12);
    EXPECT_NEAR(Pi.Integrate(ni), 1.0, 1e-12);
    const rmat_t Gi=blazem::trans(Pi.Contraction(0))*Sc*Pi.Contraction(0);
    for (size_t a=0;a<6;a++) for (size_t b=0;b<6;b++) EXPECT_NEAR(Gi(a,b), a==b?1.0:0.0, 1e-12);
    // Adjointness holds for the atomic pair as for the Löwdin one: <W, Forward(D)> == <Adjoint(W), D>.
    rvec_t W(9); for (size_t k=0;k<9;k++) W[k]=0.3-0.05*k;
    const hmat_t<double> Vw=P.Adjoint(W);
    double lhs=0; for (size_t k=0;k<9;k++) lhs+=W[k]*n[k];
    double rhs=0; for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++) rhs+=Vw(i,j)*D(j,i);
    EXPECT_NEAR(lhs, rhs, 1e-12);
}

// ORTHO-ATOMIC (slice D): two contracted manifolds on overlapping column sets are Löwdin-orthogonalised AMONG
// themselves before projecting -- the union of their functions is S-orthonormal, where the plain atomic
// projector leaves a cross-manifold overlap.  A column manifold in the same run is untouched.
TEST(LowdinProjector, OrthoAtomicManifoldsAreOrthonormalAsASet)
{
    hmat_t<double> S(6);
    for (size_t i=0;i<6;i++) S(i,i)=1.0;
    for (size_t k=0;k<3;k++) S(k,k+3)=0.4;                          // "site A" columns 0-2 overlap "site B" columns 3-5
    S(0,1)=0.1; S(3,4)=0.1;
    const std::vector<std::vector<size_t>> cols{{0,1,2},{3,4,5}};
    mat_t<double> I3(3,3,0.0); for (size_t i=0;i<3;i++) I3(i,i)=1.0;
    rmat_t Sd(6,6,0.0); for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++) Sd(i,j)=S(i,j);
    const rmat_t Sinv=blazem::inv(Sd);
    auto Wof=[&](const LowdinProjector<double>& P)                  // the union's overlap: W~ = S^{-1} T per manifold, G = W~^T S W~
    {
        rmat_t W(6,6,0.0);
        for (size_t M=0;M<2;M++) { const rmat_t w=Sinv*P.T(M); for (size_t i=0;i<6;i++) for (size_t a=0;a<3;a++) W(i,3*M+a)=w(i,a); }
        rmat_t G(6,6,0.0);
        for (size_t a=0;a<6;a++) for (size_t b=0;b<6;b++) { double t=0; for (size_t i=0;i<6;i++) for (size_t j=0;j<6;j++) t+=W(i,a)*Sd(i,j)*W(j,b); G(a,b)=t; }
        return G;
    };
    const LowdinProjector<double> atomic(S, cols, {I3,I3}, {false,false});
    const LowdinProjector<double> ortho (S, cols, {I3,I3}, {true, true});
    const rmat_t Ga=Wof(atomic), Go=Wof(ortho);
    double cross=0; for (size_t a=0;a<3;a++) for (size_t b=3;b<6;b++) cross=std::max(cross, std::abs(Ga(a,b)));
    EXPECT_GT(cross, 0.1) << "plain atomic projectors on overlapping sites are not mutually orthogonal";
    for (size_t a=0;a<6;a++) for (size_t b=0;b<6;b++) EXPECT_NEAR(Go(a,b), a==b?1.0:0.0, 1e-12);
    // The pair is still adjoint, and the coefficient map is still T^dagger.
    vec_t<double> c(6); for (size_t i=0;i<6;i++) c[i]=0.2-0.05*i;
    const vec_t<double> l1=blazem::trans(ortho.T(1))*c, l2=ortho.Coefficients(1)*c;
    for (size_t a=0;a<3;a++) EXPECT_NEAR(l1[a], l2[a], 1e-13);
    // A density that is one ortho-atomic function of manifold 1 has occupation 1 there and 0 on manifold 0.
    const vec_t<double> w=Sinv*blazem::column(ortho.T(1),0);
    hmat_t<double> D(6); for (size_t i=0;i<6;i++) for (size_t j=i;j<6;j++) D(i,j)=w[i]*w[j];
    const rvec_t nn=ortho.Forward(D);
    EXPECT_NEAR(nn[9+0], 1.0, 1e-10);                              // manifold 1 (offset 9), function 0
    double n0=0; for (size_t k=0;k<9;k++) n0+=std::abs(nn[k]);
    EXPECT_NEAR(n0, 0.0, 1e-10) << "orthogonal to every function of manifold 0";
}
