// File: Common/tests/RunPolicy.C  The CP2K-deviation policy (doc/OpenWork.md N5).
//
// What is worth gating here is the RESOLUTION RULE, not the individual defaults: the defaults are
// recorded facts that will change as routes come and go, while the rule -- explicit knob beats the
// umbrella, and a run that is not actually at parity says so -- is the thing that makes the switch
// trustworthy.  A silently-overruled knob would turn CP2K_COMPAT into a trap, which is exactly the
// failure a one-switch design invites.
#include <gtest/gtest.h>
#include <string>
#include <thread>

import qchem.RunPolicy;

using namespace qchem;

// D-ENV step 6b: the policy is the typed \c RunPolicySpec (the deck's `policy` block); no environment variable sets it any more.  Since
// D-MUTABLE-STATICS-RUNPOLICY (2026-11-24) there is no process-wide policy either: each test resolves its spec into a local value, as a
// run's setup does.

// The DEFAULT run is NOT at parity, and it is able to say which routes make it so.  This is the whole
// point of the item: the low-rank rho route has been on by default since 07d13bf6 and no timing row
// since has said so, because nothing could.
TEST(RunPolicy, DefaultRunDeviatesAndNamesTheRoutes)
{
    const RunPolicy p{RunPolicySpec{}};
    EXPECT_FALSE(p.CP2KCompat());
    EXPECT_FALSE(p.AtParity()) << "the default build runs accelerations CP2K does not -- if this ever "
                                  "passes, either a route was retired or one stopped being declared";
    EXPECT_TRUE (p.DMLowRank());
    EXPECT_TRUE (p.StreamFold());
    EXPECT_TRUE (p.MixRhoM()) << "N3 promoted 2026-09-20: (rho,m) is the default; (up,dn) is CP2K's";
    EXPECT_FALSE(p.XCFromDM());
    EXPECT_TRUE (p.SymmetryImposition()) << "by default the caller's imposeSymmetry is obeyed";
    // The collocation box tolerance: eps/|c_ij| by default, flat eps (CP2K's rule) under the umbrella.
    // It selects a LatticeScreener object, not a branch in the box walk -- doc/OldPlans/ScreeningPlan.md §7.
    EXPECT_TRUE (p.DAwareScreen());
    // DFT+U on the Löwdin block's eigenvalues (Dudarev) by default; CP2K's form is its diagonal populations.
    EXPECT_TRUE (p.HubbardEigen());
    EXPECT_EQ(p.Deviations().size(), 8u) << "a new accelerator is not finished until it is in this list";
    EXPECT_NE(p.Banner().find("DEVIATING"), std::string::npos);
}

// THE UMBRELLA: one switch turns every route off, and the policy then reports parity.
TEST(RunPolicy, CP2KCompatTurnsEveryRouteOff)
{
    RunPolicySpec spec; spec.cp2kCompat=true;
    const RunPolicy p(spec);
    EXPECT_TRUE (p.CP2KCompat());
    EXPECT_TRUE (p.AtParity());
    EXPECT_FALSE(p.DMLowRank());
    EXPECT_FALSE(p.StreamFold());
    // THE ONE THAT OVERRULES THE CALLER: CP2K does no symmetry work at all, and every banked recipe asks
    // for an imposition -- so the veto has to come from the switch, not from editing each recipe.
    EXPECT_FALSE(p.SymmetryImposition());
    // CP2K screens the collocation on geometry alone (task_list_methods.F: the radius takes no density
    // argument), so parity means the geometry-only screener.  This tree had been taking the deviation
    // SILENTLY until 2026-09-04, which is the reason it is on the table at all.
    EXPECT_FALSE(p.DAwareScreen());
    // CP2K's +U keeps only the diagonal of the Löwdin block (dft_plus_u.F) -- no eigen-decomposition.
    EXPECT_FALSE(p.HubbardEigen());
    for (const Deviation& dev : p.Deviations()) EXPECT_FALSE(dev.Deviates()) << dev.knob;
    EXPECT_NE(p.Banner().find("AT PARITY"), std::string::npos);
}

// THE OVERRIDE RULE, and the reason it matters: "parity everywhere except this one" is a legitimate
// experiment, so an explicitly-named knob outranks the umbrella -- and the run must still declare that
// it is NOT at parity, or the switch quietly lies about what was measured.
TEST(RunPolicy, AnExplicitKnobOutranksTheUmbrellaAndSaysSo)
{
    RunPolicySpec spec; spec.cp2kCompat=true; spec.streamFold=true;
    const RunPolicy p(spec);
    EXPECT_TRUE (p.CP2KCompat());
    EXPECT_TRUE (p.StreamFold())  << "an explicitly-stated knob must not be silently overruled";
    EXPECT_FALSE(p.DMLowRank())   << "...while the umbrella still governs the knobs nobody named";
    EXPECT_FALSE(p.AtParity())    << "asking for compat and then contradicting it is NOT parity";
    EXPECT_NE(p.Banner().find("(stated)"), std::string::npos);
}

// A route stated FALSE is stated.  Distinguishing "off" from "not mentioned" (an unset optional) is what makes the override rule
// expressible at all -- without it, cp2kCompat could not tell a deliberate off from a default off.
TEST(RunPolicy, StatedFalseCountsAsStated)
{
    RunPolicySpec spec; spec.streamFold=false;
    const RunPolicy p(spec);
    EXPECT_FALSE(p.StreamFold());
    bool found=false;
    for (const Deviation& d : p.Deviations())
        if (std::string(d.knob)=="policy.streamFold") {found=true; EXPECT_TRUE(d.stated); EXPECT_FALSE(d.Deviates());}
    EXPECT_TRUE(found);
    EXPECT_EQ(p.Spec(), spec) << "the policy remembers what was STATED (the record's `policy` block)";
}

// ---- D-THREADS: the ONE thread-count rule (qchem.Parallel) ----
#include "gtest/gtest.h"
import qchem.Parallel;
TEST(Parallel, ParseThreadCountIsTheOneRule)
{
    EXPECT_EQ(qchem::ParseThreadCount(nullptr, 1, 8), 1) << "unset => the default (serial)";
    EXPECT_EQ(qchem::ParseThreadCount("",      1, 8), 1) << "empty => the default";
    EXPECT_EQ(qchem::ParseThreadCount("6",     1, 8), 6);
    EXPECT_EQ(qchem::ParseThreadCount("1",     4, 8), 1) << "an explicit 1 means serial, whatever the default";
    EXPECT_EQ(qchem::ParseThreadCount("0",     1, 8), 8) << "0 means AUTO (= physical cores), in every region";
    EXPECT_EQ(qchem::ParseThreadCount("-3",    4, 8), 1);
    EXPECT_EQ(qchem::ParseThreadCount("abc",   4, 8), 1);
    EXPECT_EQ(qchem::ParseThreadCount("4x",    4, 8), 1) << "trailing garbage is not a number";
}
TEST(Parallel, PhysicalCoresIsPositiveAndNoMoreThanTheHardwareThreads)
{
    EXPECT_GE(qchem::PhysicalCores(), 1);
    EXPECT_LE(unsigned(qchem::PhysicalCores()), std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 1u);
}
TEST(Parallel, ThreadSummaryNamesTheKnobs)
{
    const std::string s=qchem::ThreadSummary();
    EXPECT_NE(s.find("OpenMP regions="), std::string::npos);
    EXPECT_NE(s.find("BLAS="),           std::string::npos);
}
