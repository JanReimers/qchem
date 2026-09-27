// File: Response/tests/Ring.C  chi0 of an exactly solvable model, against brute force (stage R0).
//
// THE MODEL.  A ring of N cells, two sites per cell (a, b), one orbital per site, spinless:
//     <a,R|H|a,R> = ea,  <b,R|H|b,R> = eb,  <a,R|H|b,R> = t,  <a,R+1|H|b,R> = t'  (and Hermitian partners).
// Its N Bloch blocks, in the basis gauge the library uses (phase e^{ik.R} over lattice translations only),
// are the 2x2 matrices H_ab(k) = t + t' e^{-2 pi i k}.  The channels are the two SITES: P_J = |J><J|, whose
// amplitudes on block k are just row J of the eigenvector matrix.
//
// THE ORACLE.  The ring's own 2N x 2N Hamiltonian, diagonalised directly, perturbed by alpha |J,0><J,0| on
// ONE cell, refilled, and the occupations n_{I,R} measured by central difference -- the non-self-consistent
// response of an ISOLATED perturbation.  The N-point q-mesh on the N-point k-mesh is mathematically that
// ring, so ChannelResponse::RealSpace() must reproduce every element, off-diagonal and inter-cell ones
// included: this pins the q-pairing, the phase convention, the Fourier sum and (for the metal) the q=0
// Fermi shift δμ that keeps the electron count fixed.  No SCF, no basis, no Hamiltonian library in the loop.
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <vector>
import qchem.Response.Probe;
import qchem.Symmetry.Factory;   // BlochFactory
import qchem.Symmetry.Irrep;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::Response;

namespace {

double Fermi(double e, double mu, double kT) {return 1.0/(1.0+std::exp((e-mu)/kT));}

//! The chemical potential that puts \a nel electrons in levels \a e (bisection; kT>0).
double SolveMu(const std::vector<double>& e, double nel, double kT)
{
    double lo=-50, hi=50;
    for (int it=0;it<200;it++)
    {
        const double mu=0.5*(lo+hi);
        double n=0; for (double x : e) n+=Fermi(x,mu,kT);
        (n>nel ? hi : lo)=mu;
    }
    return 0.5*(lo+hi);
}

struct Ring
{
    int    N=6;
    double ea=-0.3, eb=0.4, t=0.5, tp=0.25;
    double kT=0.0;          //!< 0 = integer (insulator: one electron per cell, per-block fill)
    double nel=1.0;         //!< electrons per cell (a metal takes a fraction under Fermi smearing)

    //! Occupations of the brute-force ring's levels (sorted ascending), total N*nel electrons.
    std::vector<double> Fill(const rvec_t& e) const
    {
        std::vector<double> f(e.size(),0.0);
        if (kT==0.0) { for (int i=0;i<N*int(nel);i++) f[i]=1.0; return f; }
        std::vector<double> ev; for (size_t i=0;i<e.size();i++) ev.push_back(e[i]);   // no std range ctor on blaze iterators
        const double mu=SolveMu(ev, N*nel, kT);
        for (size_t i=0;i<e.size();i++) f[i]=Fermi(e[i],mu,kT);
        return f;
    }

    //! n_{I,R} of the ring perturbed by alpha on site J of cell 0.
    std::vector<double> Occupations(int J, double alpha) const
    {
        const int M=2*N;
        rsmat_t H(M);
        for (int R=0;R<N;R++)
        {
            H(2*R,2*R)=ea; H(2*R+1,2*R+1)=eb;
            H(2*R,2*R+1)=t;
            H(2*((R+1)%N),2*R+1)=tp;     // <a,R+1|H|b,R>
        }
        H(J,J)+=alpha;                   // site J of cell 0
        rvec_t e; rmat_t U;
        blazem::eigen(H,e,U);
        const std::vector<double> f=Fill(e);
        std::vector<double> n(M,0.0);
        for (int i=0;i<M;i++) for (int s=0;s<M;s++) n[s]+=f[i]*U(s,i)*U(s,i);
        return n;
    }

    //! The Bloch reference + the two site channels (amplitude = row J of the eigenvectors).
    struct Built { std::unique_ptr<Reference> ref; std::vector<std::vector<cmat_t>> amp; };
    Built Build() const
    {
        std::vector<ReferenceBlock> blocks;
        std::vector<std::vector<cmat_t>> amp;
        std::vector<double> allE;
        std::vector<rvec_t> es;
        for (int ik=0;ik<N;ik++)
        {
            const double k=double(ik)/N;
            chmat_t H(2);
            H(0,0)=ea; H(1,1)=eb;
            H(0,1)=t+tp*std::exp(dcmplx(0.0,-2.0*M_PI*k));
            rvec_t e; cmat_t C;
            blazem::eigen(H,e,C);
            ReferenceBlock b;
            b.irrep=Irrep(Spin::Up, Symmetry::BlochFactory(ivec3_t(N,1,1),ivec3_t(ik,0,0),1.0/N));
            b.w=1.0/N; b.g=1.0; b.e=e; b.f=rvec_t(2,0.0); b.reservoir=0;
            blocks.push_back(b);
            es.push_back(e);
            for (size_t i=0;i<e.size();i++) allE.push_back(e[i]);
            std::vector<cmat_t> a;
            for (int J=0;J<2;J++) { cmat_t L(1,2); L(0,0)=C(J,0); L(0,1)=C(J,1); a.push_back(L); }
            amp.push_back(a);
        }
        if (kT==0.0) for (auto& b : blocks) {b.f[0]=1.0; b.f[1]=0.0;}   // one electron per block (per-block fill)
        else
        {
            const double mu=SolveMu(allE, N*nel, kT);
            for (auto& b : blocks) for (size_t n=0;n<2;n++) b.f[n]=Fermi(b.e[n],mu,kT);
        }
        OccupationConfig cfg; cfg.kT=kT;
        return {std::make_unique<Reference>(std::move(blocks), MakeOccupancyRule(cfg), 1e-12), std::move(amp)};
    }

    void CheckAgainstBruteForce(double tol) const
    {
        Built B=Build();
        AmplitudeProbe probe(*B.ref, B.amp, {"a","b"});
        auto chi=IndependentResponse(*B.ref, probe, ivec3_t(N,1,1));
        ASSERT_TRUE(chi.IsOk()) << (chi ? "" : chi.Error().detail);
        const rmat_t R=chi->RealSpace();
        const double h=1e-5;
        for (int J=0;J<2;J++)
        {
            const std::vector<double> np=Occupations(J,+h), nm=Occupations(J,-h);
            for (int r=0;r<N;r++)
                for (int I=0;I<2;I++)
                {
                    const double fd=(np[2*r+I]-nm[2*r+I])/(2*h);
                    EXPECT_NEAR(R(r*2+I, J), fd, tol) << "I=" << I << " R=" << r << " J=" << J;
                }
        }
    }
};

} // namespace

TEST(ResponseRing, InsulatorChi0MatchesBruteForceRing)
{
    Ring().CheckAgainstBruteForce(1e-7);
}

TEST(ResponseRing, MetalChi0WithFermiShiftMatchesBruteForceRing)
{
    Ring m;
    m.ea=0.0; m.eb=0.1; m.t=0.4; m.tp=0.3; m.kT=0.05; m.nel=0.7; m.N=8;
    m.CheckAgainstBruteForce(1e-7);
}

TEST(ResponseRing, ChannelResponseIsHermitianAndRealSpaceIsSymmetric)
{
    Ring ring;
    auto B=ring.Build();
    AmplitudeProbe probe(*B.ref, B.amp, {"a","b"});
    auto chi=IndependentResponse(*B.ref, probe, ivec3_t(ring.N,1,1));
    ASSERT_TRUE(chi.IsOk());
    for (const auto& c : chi->chi)
        for (size_t I=0;I<2;I++) for (size_t J=0;J<2;J++)
            EXPECT_NEAR(std::abs(c(I,J)-std::conj(c(J,I))), 0.0, 1e-12);
    const rmat_t R=chi->RealSpace();
    for (size_t i=0;i<R.rows();i++) for (size_t j=0;j<R.columns();j++) EXPECT_NEAR(R(i,j), R(j,i), 1e-12);
    // An insulator's electron count cannot respond: every column of the real-space response sums to zero.
    for (size_t j=0;j<R.columns();j++) { double s=0; for (size_t i=0;i<R.rows();i++) s+=R(i,j); EXPECT_NEAR(s,0.0,1e-12); }
    EXPECT_LT(R(0,0), 0.0);   // a potential raised on a site drains it
}

TEST(ResponseRing, IncommensurateQMeshFails)
{
    Ring ring;
    auto B=ring.Build();
    AmplitudeProbe probe(*B.ref, B.amp, {"a","b"});
    auto chi=IndependentResponse(*B.ref, probe, ivec3_t(4,1,1));   // 4 does not divide 6
    ASSERT_FALSE(chi.IsOk());
    EXPECT_EQ(chi.Error().why, ResponseFailure::Why::Incommensurate);
}

// E1's review contract: an integer per-block fill that is NOT aufbau across a coupled (k, k+q) pair is a
// FAILED Outcome naming the pair -- never a huge weight that silently corrupts chi.  Two k-blocks, each
// with one electron: the occupied level at k=1/2 (0.5) sits ABOVE the empty level at k=0 (0.2).
namespace {
std::unique_ptr<Reference> TwoBlock(double eOcc1, double eEmp0, double noise)
{
    std::vector<ReferenceBlock> blocks(2);
    for (int ik=0;ik<2;ik++)
    {
        blocks[ik].irrep=Irrep(Spin::Up, Symmetry::BlochFactory(ivec3_t(2,1,1),ivec3_t(ik,0,0),0.5));
        blocks[ik].w=0.5; blocks[ik].g=1.0; blocks[ik].f={1.0,0.0};
    }
    blocks[0].e={-0.4, eEmp0};
    blocks[1].e={eOcc1, 0.9};
    return std::make_unique<Reference>(std::move(blocks), MakeOccupancyRule(OccupationConfig{}), noise);
}
AmplitudeProbe OneChannel(const Reference& ref)
{
    std::vector<std::vector<cmat_t>> amp(2);
    for (auto& a : amp) { cmat_t L(1,2); L(0,0)=0.8; L(0,1)=0.6; a.push_back(L); }
    return AmplitudeProbe(ref, amp, {"x"});
}
}

TEST(ResponseRing, InvertedCoupledPairIsAFailedOutcome)
{
    auto ref=TwoBlock(0.5, 0.2, 1e-8);
    auto chi=IndependentResponse(*ref, OneChannel(*ref), ivec3_t(2,1,1));
    ASSERT_FALSE(chi.IsOk());
    EXPECT_EQ(chi.Error().why, ResponseFailure::Why::Inverted);
    EXPECT_NE(chi.Error().detail.find("INVERTED"), std::string::npos);
    // q = 0 alone couples each block with itself, where the per-block fill IS aufbau: fine.
    auto q0=IndependentResponse(*ref, OneChannel(*ref), ivec3_t(1,1,1));
    EXPECT_TRUE(q0.IsOk());
}

TEST(ResponseRing, UnresolvedGapIsAFailedOutcomeButASmallResolvedGapIsPhysics)
{
    // Gap 5e-4 Ha across the q pair, eigenvalue noise 1e-3: unresolved -> fail.
    auto bad=TwoBlock(0.1995, 0.2, 1e-3);
    auto c1=IndependentResponse(*bad, OneChannel(*bad), ivec3_t(2,1,1));
    ASSERT_FALSE(c1.IsOk());
    EXPECT_EQ(c1.Error().why, ResponseFailure::Why::Unresolved);
    // The SAME small gap with noise 1e-7 is resolved: a large chi0, reported with its bound, not an error.
    auto ok=TwoBlock(0.1995, 0.2, 1e-7);
    auto c2=IndependentResponse(*ok, OneChannel(*ok), ivec3_t(2,1,1));
    ASSERT_TRUE(c2.IsOk());
    EXPECT_NEAR(c2->gap, 5e-4, 1e-12);
}

// A recipe that computes no [F,D] hands the reference an UNMEASURED noise (NaN).  The gate must still catch
// an inverted pair (the sign needs no noise), must NOT invent an "unresolved" verdict, and must say
// "UNMEASURED" rather than print a false 0.
TEST(ResponseRing, UnmeasuredNoiseGatesTheSignOnlyAndSaysSo)
{
    const double nan=std::numeric_limits<double>::quiet_NaN();
    auto inv=TwoBlock(0.5, 0.2, nan);
    auto c1=IndependentResponse(*inv, OneChannel(*inv), ivec3_t(2,1,1));
    ASSERT_FALSE(c1.IsOk());
    EXPECT_EQ(c1.Error().why, ResponseFailure::Why::Inverted);
    EXPECT_NE(c1.Error().detail.find("UNMEASURED"), std::string::npos);
    auto tiny=TwoBlock(0.1995, 0.2, nan);                  // a 5e-4 gap: resolved or not is unknowable
    auto c2=IndependentResponse(*tiny, OneChannel(*tiny), ivec3_t(2,1,1));
    EXPECT_TRUE(c2.IsOk());
}
