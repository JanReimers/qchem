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

// A Fock with a clear 2-occupied / 2-virtual split, plus off-diagonal occ-virt coupling so the orbital
// gradient is NONZERO -- otherwise ComputeStep has nothing to step along and the test is vacuous.
rsmat_t MakeFock()
{
    rsmat_t F(N);
    F(0,0)=-2.0; F(1,1)=-1.5; F(2,2)=1.0; F(3,3)=1.5;
    F(0,2)=0.20; F(1,3)=0.15;      // occ-virt coupling => nonzero gradient
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

// A GDM accelerator seeded far enough to hold a live geodesic step: one diagonalizing NextOrbitals (which
// caches the orbitals), then a second UseFD so the residual is current.  FDMax is set high so the [F,D]
// gate never blocks the step -- this test is about rejection, not about engagement.
struct Rig
{
    GDMParams                    params{ .FDMax=1e30, .Trust=0.1, .TrustBackoff=0.25, .TrustMin=1e-4 };
    LASolver<double>*            las = LASolver<double>::Factory(qchem::Cholesky);
    tSCFAcceleratorGDM<double>   acc{params};
    tSCFIrrepAccelerator<double>* irrep=nullptr;

    Rig()
    {
        las->SetBasisOverlap(Identity(N));
        irrep = acc.Create(las, Irrep(), (int)NOCC);
        const rsmat_t F=MakeFock();
        irrep->UseFD(F, MakeDPrime(rmat_t(Identity(N))));
        auto [U,Up,e] = irrep->NextOrbitals();     // diagonalize: caches orbitals (itsHaveC)
        irrep->UseFD(F, MakeDPrime(Up));           // fresh residual at the new orbitals
    }
    ~Rig() { delete las; }   // acc owns the irrep accelerators
};
} //anon

// The rig itself must hold a LIVE step, or every assertion below passes vacuously.
TEST(SCFAcceleratorGDM, RigHoldsALiveGeodesicStep)
{
    Rig r;
    EXPECT_TRUE(r.acc.CanLineSearch()) << "the [F,D] gate should be open (FDMax=1e30)";
    EXPECT_TRUE(r.irrep->ComputeStep()) << "a seeded GDM with a nonzero gradient must offer a step";
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
    r.irrep->UseFD(MakeFock(), MakeDPrime(Up));
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
    r.irrep->UseFD(MakeFock(), MakeDPrime(Up0));
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
