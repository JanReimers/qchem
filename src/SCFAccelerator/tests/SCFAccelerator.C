// File: SCFAccelerator/tests/SCFAccelerator.C  Unit tests for the accelerator STEP-REJECTION contract.
//
// An accelerator PROPOSES a step; only the caller can judge it, because only the caller can evaluate the
// energy (the accelerator sits below the wavefunction/Hamiltonian in the library DAG).  RejectStep() is how
// that judgement comes back.  The contract has one clause that is a SAFETY guarantee rather than a quality
// knob, and it is the clause worth a test:
//
//     RejectStep()==false ("exhausted") MUST leave the accelerator in a state where ComputeStep() is false
//     and CanLineSearch() is false.
//
// Without it the caller cannot fall back at all: its fallback runs NextOrbitals(), and a minimizer still
// holding a live step TAKES that step there (GDM: OrbitalsAt(1.0,true)) -- the very step just rejected, at
// FULL length.  A regression in that clause would not throw or crash; it would silently re-commit rejected
// steps, which is exactly the defect these tests were written for (doc/SymmetryUpgradePlan.md §7 step 7:
// a committed non-descent step threw MnO by +14.5 Ha).
//
// Mock setup -- no BasisSet, no Hamiltonian: an ORTHONORMAL basis (S=I) so the ortho transform is the
// identity, and a diagonally-dominant Fock with a clean occ/virt gap.
#include <gtest/gtest.h>
import qchem.SCFAccelerator;
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
import qchem.LASolver;
import qchem.Symmetry;
import qchem.Blaze;
using namespace qchem;
using namespace qchem::SCFAccelerators;

namespace
{
constexpr size_t N=4, NOCC=2;

rsmat_t Identity(size_t n) { rsmat_t I(n); for (size_t i=0;i<n;i++) I(i,i)=1.0; return I; }

// TWO Focks, and the second one is the whole point.  A rig that diagonalizes F and then asks GDM to step
// on that SAME F is at the exact minimum already: in its own eigenbasis the occupied-virtual block of F is
// zero, so the gradient vanishes and every step is a no-op.  (Measured: the band energy was bit-identical
// at t=0, 0.05 and 1.0 -- which is what exposed the first draft of these tests as vacuous.)  A real SCF
// iteration never looks like that: the orbitals come from the PREVIOUS Fock and the current one has moved.
// So: seed by diagonalizing F1, then hand the accelerator F2.  The o-v gradient is then genuinely nonzero
// and there is something to descend.
rsmat_t MakeFock()
{
    rsmat_t F(N);
    F(0,0)=-2.0; F(1,1)=-1.5; F(2,2)=1.0; F(3,3)=1.5;
    F(0,2)=0.20; F(1,3)=0.15;      // occ-virt coupling => nonzero gradient
    return F;
}

// The SECOND Fock: same gapped structure, different couplings, so F1's eigenvectors are not F2's.
rsmat_t MakeFock2()
{
    rsmat_t F(N);
    F(0,0)=-2.1; F(1,1)=-1.4; F(2,2)=1.1; F(3,3)=1.6;
    F(0,2)=0.35; F(1,3)=0.05; F(0,3)=0.10;
    return F;
}

// D' for the lowest NOCC orbitals of the CURRENT orbital set (columns of C).
rsmat_t MakeDPrime(const rmat_t& C)
{
    rsmat_t D(N);
    for (size_t i=0;i<N;i++)
        for (size_t j=0;j<=i;j++)
        {
            double v=0.0;
            for (size_t k=0;k<NOCC;k++) v += C(i,k)*C(j,k);
            D(i,j)=v;
        }
    return D;
}

// The band energy of the occupied block, E(C)=sum_k <c_k|F|c_k> over the leading NOCC columns -- the
// quantity a direct minimizer descends at FIXED F, so a correct search direction must lower it for
// small enough t.  That is precisely the property a sign error in the preconditioner destroys: PG
// divides the gradient by (eps_a - eps_i), and an unfloored gap flips sign on a non-aufbau pair,
// turning part of "steepest descent" into ascent (hence GDMParams::PCFloor).  The rest of this file
// tests the bail-out BOOKKEEPING; this is what tests that the thing being bailed out of points the
// right way.
double BandEnergy(const rmat_t& Cp, const rsmat_t& F)
{
    double e=0.0;
    for (size_t k=0;k<NOCC;k++)
        for (size_t i=0;i<N;i++)
            for (size_t j=0;j<N;j++)
                e += Cp(i,k)*F(i,j)*Cp(j,k);
    return e;
}

// A GDM accelerator seeded far enough to hold a live geodesic step: one diagonalizing NextOrbitals (which
// caches the orbitals), then a second UseFD so the residual is current.  FDMax is set high so the [F,D]
// gate never blocks the step -- this test is about rejection, not about engagement.
struct Rig
{
    GDMParams                    params{ .FDMax=1e30, .Trust=0.1, .TrustBackoff=0.25, .TrustMin=1e-4 };
    LASolver<double>*            las = LASolver<double>::Factory(qchem::Cholesky);
    tSCFAcceleratorGDM<double>   acc{params};
    tSCFIrrepAccelerator<double>* irrep=nullptr;

    rsmat_t F2=MakeFock2();   // the Fock the accelerator is CURRENTLY on (see MakeFock's note)
    Rig()
    {
        las->SetBasisOverlap(Identity(N));
        irrep = acc.Create(las, Irrep(), (int)NOCC);
        irrep->UseFD(MakeFock(), MakeDPrime(rmat_t(Identity(N))));
        auto [U,Up,e] = irrep->NextOrbitals();     // diagonalize F1: caches orbitals (itsHaveC)
        irrep->UseFD(F2, MakeDPrime(Up));          // ...then MOVE to F2, so the o-v gradient is nonzero
    }
    ~Rig() { delete las; }   // acc owns the irrep accelerators
};
} //anon

// The rig itself must hold a LIVE step, or every assertion below passes vacuously.
TEST(SCFAcceleratorGDM, RigHoldsALiveGeodesicStep)
{
    Rig r;
    EXPECT_TRUE(r.acc.CanLineSearch()) << "the [F,D] gate should be open (FDMax=1e30)";
    EXPECT_TRUE(r.irrep->ComputeStep()) << "a seeded GDM must offer a step";
    // ...and the step must be REAL: on a Fock the orbitals do not diagonalize, the band energy
    // must actually move.  Without this the whole file could pass on a zero gradient.
    auto [U0,Cp0,e0] = r.irrep->OrbitalsAt(0.0,false);
    auto [U1,Cp1,e1] = r.irrep->OrbitalsAt(1.0,false);
    EXPECT_GT(std::abs(BandEnergy(Cp1,r.F2)-BandEnergy(Cp0,r.F2)), 1e-12)
        << "the rig is at a stationary point -- every other test here would be vacuous";
}

// THE SAFETY CLAUSE.  Reject until exhausted, then check the guarantee the caller relies on.
TEST(SCFAcceleratorGDM, ExhaustedRejectionForcesADiagonalizeSoTheCallerCanFallBackSafely)
{
    Rig r;
    ASSERT_TRUE(r.irrep->ComputeStep());
    // Trust 0.1 backing off by 0.25 hits TrustMin=1e-4 after 5 rejections; allow slack but require
    // TERMINATION -- an accelerator that never exhausts would spin the caller's retry loop forever.
    int rejects=0;
    while (r.acc.RejectStep()) { ASSERT_LT(++rejects, 50) << "RejectStep never reports exhaustion"; }
    EXPECT_GT(rejects, 0) << "the trust radius should buy at least one retry before giving up";

    EXPECT_FALSE(r.irrep->ComputeStep())
        << "EXHAUSTED must silence ComputeStep -- else the caller's NextOrbitals() fallback re-takes the "
           "rejected step at full length";
    EXPECT_FALSE(r.acc.CanLineSearch())
        << "EXHAUSTED must also close the line-search gate, so the iterator picks the fixed-point driver";
}

// ...and the forced diagonalize must be CONSUMED by an actual diagonalize, not by a mere query -- otherwise
// a CanLineSearch()/ComputeStep() probe (the iterator makes several per iteration) would clear the guarantee
// before the fallback step it was armed for ever runs.
TEST(SCFAcceleratorGDM, TheForcedDiagonalizeSurvivesQueriesAndIsClearedOnlyByADiagonalize)
{
    Rig r;
    ASSERT_TRUE(r.irrep->ComputeStep());
    while (r.acc.RejectStep()) {}

    for (int i=0;i<3;i++)                                  // queries must NOT consume the flag
    {
        EXPECT_FALSE(r.irrep->ComputeStep()) << "query " << i << " consumed the forced diagonalize";
        EXPECT_FALSE(r.acc.CanLineSearch())  << "query " << i << " consumed the forced diagonalize";
    }

    auto [U,Up,e] = r.irrep->NextOrbitals();               // THE diagonalize the flag was armed for
    r.irrep->UseFD(r.F2, MakeDPrime(Up));
    EXPECT_TRUE(r.acc.CanLineSearch()) << "after diagonalizing, the geodesic must be re-armed";
    EXPECT_TRUE(r.irrep->ComputeStep()) << "a rejection must not disable the minimizer permanently";
}

// A rejection is a per-step backoff, not a permanent ratchet: an ACCEPTED step restores the full trust
// radius, so one bad direction cannot leave the minimizer crawling for the rest of the run.  Probed
// behaviourally (the trust radius is private): after rejections shrink it, a committed step must restore
// the ability to absorb many further rejections -- i.e. the backoff budget is replenished.
TEST(SCFAcceleratorGDM, AnAcceptedStepRestoresTheTrustRadius)
{
    Rig r;
    ASSERT_TRUE(r.irrep->ComputeStep());
    int firstBudget=0;
    while (r.acc.RejectStep()) ++firstBudget;
    ASSERT_GT(firstBudget, 0);

    auto [U0,Up0,e0] = r.irrep->NextOrbitals();            // clear the forced diagonalize
    r.irrep->UseFD(r.F2, MakeDPrime(Up0));
    ASSERT_TRUE(r.irrep->ComputeStep());
    r.irrep->OrbitalsAt(1.0,/*commit*/true);               // an ACCEPTED step
    ASSERT_TRUE(r.irrep->ComputeStep());

    int secondBudget=0;
    while (r.acc.RejectStep()) ++secondBudget;
    EXPECT_EQ(secondBudget, firstBudget)
        << "an accepted step must restore the full trust radius (backoff is recovery, not a ratchet)";
}

// Every non-minimizer inherits the base default: nothing to reject, no retry to offer.  That default is what
// makes the caller's fallback trivially safe for DIIS/Null (their ComputeStep() is false to begin with).
TEST(SCFAccelerator, NonMinimizersDeclineToRetryByDefault)
{
    struct Plain : public virtual tSCFIrrepAccelerator<double>
    {
        void UseFD(const hmat_t<double>&, const hmat_t<double>&) override {}
        LASolver<double>::UUd_t NextOrbitals() override { return {}; }
    } p;
    EXPECT_FALSE(p.ComputeStep()) << "a non-minimizer holds no step";
    EXPECT_FALSE(p.RejectStep())  << "...so it can offer no retry";
}


TEST(SCFAcceleratorGDM, TheGeodesicStepDescendsTheBandEnergy)
{
    Rig r;
    ASSERT_TRUE(r.irrep->ComputeStep());

    auto [U0,Cp0,e0] = r.irrep->OrbitalsAt(0.0,/*commit*/false);   // t=0 IS the pre-step point
    const double E0=BandEnergy(Cp0,r.F2);
    // Small t: a geodesic is a rotation, so a descent direction is only guaranteed to descend for
    // sufficiently small steps -- test the guarantee, not the full model step.
    auto [Ut,Cpt,et] = r.irrep->OrbitalsAt(0.05,/*commit*/false);
    EXPECT_LT(BandEnergy(Cpt,r.F2), E0)
        << "a geodesic step must LOWER the occupied-block band energy at fixed F";
}

// ...and the floor must not have broken the ordinary well-gapped case by over-damping it: the model step
// (t=1) should still descend when the occ/virt gap is healthy, the regime PCFloor=0.1 sits below.
TEST(SCFAcceleratorGDM, TheModelStepStillDescendsOnAWellGappedProblem)
{
    Rig r;
    ASSERT_TRUE(r.irrep->ComputeStep());
    auto [U0,Cp0,e0] = r.irrep->OrbitalsAt(0.0,false);
    auto [U1,Cp1,e1] = r.irrep->OrbitalsAt(1.0,false);
    EXPECT_LT(BandEnergy(Cp1,r.F2), BandEnergy(Cp0,r.F2))
        << "PCFloor must not damp a healthy gap (this Fock's occ/virt gap is ~1 Ha, PCFloor is 0.1)";
}
