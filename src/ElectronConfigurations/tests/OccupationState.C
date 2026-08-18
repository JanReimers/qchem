// File: ElectronConfigurations/tests/OccupationState.C  The occupation seam's persistent memory (R2.21).
//
// What these pin is the SPLIT itself, not the fill arithmetic (that is exercised end-to-end by the SCF
// suite).  Three properties, each one a defect the split repairs:
//
//  1. ONE LEDGER, MANY SCALARS.  A mixed real/complex mesh runs one entropy aggregate, one fill clock and
//     one arming flag across blocks of BOTH scalars.  Before the split the ledger lived inside a
//     run-typed policy, so a <double> block in a <dcmplx> run could only reach it through a forwarding
//     view -- and forwarding is the kind of thing that is correct until someone adds a field.
//  2. REFERENCES ARE TYPED PER BLOCK.  A real block's MOM reference is a mat_t<double>, a complex block's
//     a mat_t<dcmplx>, and both live in the same map.  This is the one that unblocked MOM on a mixed
//     mesh: a run-typed reference home had nowhere to put the real block's matrix, so its capture threw.
//  3. A SCALAR MIX-UP IS LOUD.  Reading a block's reference at the wrong scalar THROWS rather than
//     reading as "no reference yet" -- which would silently restart the delayed-IMOM capture clock.
#include "gtest/gtest.h"
#include <complex>
#include <vector>
import qchem.ElectronConfiguration.OccupationPolicy;
import qchem.Symmetry.Factory;   // YFactory (a concrete Irrep to key blocks by)
import qchem.Symmetry.Irrep;
import qchem.Types;

using namespace qchem;

namespace {

//! A minimal OrbitalView: \a nocc majority-filled orbitals of length \a n, each a distinct unit vector,
//! so an adopted reference is an identity-like column set whose overlaps are exactly known.
template <class T> class FakeOrbitals : public OrbitalView<T>
{
public:
    FakeOrbitals(size_t n, size_t nocc) : itsN(n), itsNocc(nocc) {}
    size_t   NumOrbitals()        const override {return itsN;}
    double   Occupation (size_t i) const override {return i<itsNocc ? 1.0 : 0.0;}
    double   Degeneracy (size_t)   const override {return 1.0;}
    vec_t<T> CoeffPrime (size_t i) const override
    { vec_t<T> c(itsN, T(0.0)); c[i]=T(1.0); return c; }
private:
    size_t itsN, itsNocc;
};

Irrep S(Spin ms) {return Irrep(ms, Symmetry::YFactory(0));}
Irrep P(Spin ms) {return Irrep(ms, Symmetry::YFactory(1));}

} // anon

// (1) The ledger is scalar-INDEPENDENT: entropy, fill clock and arming are shared by every policy over
// one state, whatever scalar each one serves.  This is what makes a mixed mesh's bookkeeping
// indistinguishable from its all-complex twin's.
TEST(OccupationState, OneLedgerIsSharedAcrossScalars)
{
    OccupationState st;
    OccupationPolicy<dcmplx> cplx(st);
    OccupationPolicy<double> real_(st);

    cplx.BeginFill();
    cplx.AccumulateEntropy(-0.25);
    real_.AccumulateEntropy(-0.75);          // the REAL block's -TS lands in the same aggregate
    EXPECT_DOUBLE_EQ(cplx.EntropyTerm(), -1.0);
    EXPECT_DOUBLE_EQ(real_.EntropyTerm(), -1.0) << "both policies must read ONE run's -TS";

    real_.BeginFill();                        // either side may open the next fill
    EXPECT_DOUBLE_EQ(cplx.EntropyTerm(), 0.0);

    real_.CountFill(S(Spin::Up));             // one clock per BLOCK, shared across policies
    cplx.CountFill(S(Spin::Up));
    EXPECT_EQ(st.FillCount(S(Spin::Up)), 2);
    EXPECT_EQ(st.FillCount(P(Spin::Up)), 0)   << "an untouched block's clock stays at zero";

    EXPECT_FALSE(cplx.CrossIrrepMOMArmed());
    st.ArmCrossIrrepMOM();
    EXPECT_TRUE(real_.CrossIrrepMOMArmed()) << "arming is a property of the RUN, not of one scalar";
}

// (2) THE POINT OF R2.21: blocks of different scalars hold references in one ledger, each under its own
// type.  This is the mixed-mesh MOM case -- a real TRIM block beside complex general-k blocks.
TEST(OccupationState, ReferencesAreTypedPerBlockAndCoexist)
{
    OccupationState st;
    const Irrep gamma=S(Spin::Up), kpt=P(Spin::Up);

    st.AdoptReference(gamma, FakeOrbitals<double>(4,2));   // the real TRIM block
    st.AdoptReference(kpt,   FakeOrbitals<dcmplx>(4,3));   // a complex general-k block

    EXPECT_TRUE(st.HasReference(gamma));
    EXPECT_TRUE(st.HasReference(kpt));
    ASSERT_NE(st.Reference<double>(gamma), nullptr);
    ASSERT_NE(st.Reference<dcmplx>(kpt),   nullptr);
    EXPECT_EQ(st.Reference<double>(gamma)->columns(), 2u);  // the majority-filled columns, that scalar
    EXPECT_EQ(st.Reference<dcmplx>(kpt)->columns(),   3u);

    // Scores read each block at its own scalar: an identity-like reference scores its own occupied
    // orbitals 1 and the virtuals 0 -- the property the MOM ranking is built on.
    rvec_t s=st.Scores(gamma, FakeOrbitals<double>(4,2));
    ASSERT_EQ(s.size(), 4u);
    EXPECT_NEAR(s[0], 1.0, 1e-14);
    EXPECT_NEAR(s[1], 1.0, 1e-14);
    EXPECT_NEAR(s[2], 0.0, 1e-14);

    st.ReleaseReferences();                                 // the 0h guard drops BOTH and restarts the clocks
    EXPECT_FALSE(st.HasReference(gamma));
    EXPECT_FALSE(st.HasReference(kpt));
    EXPECT_EQ(st.FillCount(gamma), 0);
}

// (3) A wrong-scalar read is a DEFECT, not a miss.  Returning null would read as "no reference yet" and
// silently restart the delayed-IMOM capture clock -- a wrong answer that looks like a fresh run.
TEST(OccupationState, WrongScalarReferenceThrows)
{
    OccupationState st;
    const Irrep q=S(Spin::Up);
    st.AdoptReference(q, FakeOrbitals<double>(3,1));
    EXPECT_NE(st.Reference<double>(q), nullptr);
    EXPECT_THROW((void)st.Reference<dcmplx>(q), std::logic_error);
    EXPECT_EQ(st.Reference<double>(S(Spin::Down)), nullptr) << "an untouched block is a plain miss, not a throw";
}

// A capture with nothing majority-filled holds NO reference (rather than an empty matrix that reads as
// "referenced"): under Fermi smearing every orbital carries a sliver of occupation, so the majority test
// is what keeps the MOM overlap meaningful.
TEST(OccupationState, NothingMajorityFilledHoldsNoReference)
{
    OccupationState st;
    const Irrep q=S(Spin::Up);
    st.AdoptReference(q, FakeOrbitals<double>(4,0));
    EXPECT_FALSE(st.HasReference(q));
    EXPECT_TRUE(st.Scores(q, FakeOrbitals<double>(4,0)).size()==0u);
}

// The delayed-IMOM capture is the POLICY's decision over the STATE's clock: no capture before the
// configured start iteration, exactly one after -- and none at all when MOM is off.
TEST(OccupationState, DelayedCaptureIsThePolicyDecisionOverTheSharedClock)
{
    OccupationState st;
    OccupationPolicy<double> pol(st);
    const Irrep q=S(Spin::Up);
    FakeOrbitals<double> orbs(4,2);

    pol.Configure(/*useMOM*/false, /*startIter*/2, /*kT*/0.0, /*penalty*/0.0);
    pol.CountFill(q); pol.CountFill(q); pol.CountFill(q);
    pol.CaptureReferenceIfDue(q, orbs);
    EXPECT_FALSE(st.HasReference(q)) << "MOM off: never capture";

    OccupationState st2;
    OccupationPolicy<double> pol2(st2);
    pol2.Configure(true, 2, 0.0, 0.0);
    pol2.CountFill(q);
    pol2.CaptureReferenceIfDue(q, orbs);
    EXPECT_FALSE(st2.HasReference(q)) << "before the settling delay: no reference yet";
    pol2.CountFill(q);
    pol2.CaptureReferenceIfDue(q, orbs);
    EXPECT_TRUE(st2.HasReference(q))  << "at the start iteration: captured once";
}
