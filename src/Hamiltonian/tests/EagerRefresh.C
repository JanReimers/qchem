// File: Hamiltonian/tests/EagerRefresh.C  The eager refresh phase, pinned as a contract.
//
// doc/OpenWork.md item **KP** (2026-09-07).  A term's expensive density-derived state -- the Hartree
// \f$V_H\f$ field, \f$\rho\f$ on the XC quadrature's points -- is k-INDEPENDENT: one object, correct for
// every Bloch block of an SCF iteration.  It used to be filled LAZILY by whichever block asked first,
// which makes a read-only shared resource into a WRITE-ON-FIRST-TOUCH and is the one thing standing
// between `tCompositeWF::DoSCFIteration`'s per-block loop and running its blocks concurrently.
//
// `tHamiltonian::RefreshForDensity` hoists that fill into an explicit phase before the loop.  THIS FILE
// TESTS THE PLUMBING, which is the part that can silently rot: that the Hamiltonian folds the call over
// every DYNAMIC term exactly once, that a STATIC term is never asked (it is density-independent by
// definition, so it has nothing to warm), and and -- since 2026-09-10 -- that EVERY term answers, because both hooks are pure (a defaulted no-op
// was letting terms that needed the phase silently skip it).
//
// ⚠ WHY A SPY AND NOT A REAL RUN.  The production terms warm real caches whose freshness is only
// observable through a full GPW SCF; asserting "the second call was a cache hit" there needs an
// integration test and a timing argument.  The fold itself is pure plumbing over the term lists, and
// plumbing is exactly what a unit test can pin exactly -- CALL COUNTS, not timings.  (User rule: the dev
// loop is unit, integration is acceptance only.)
#include "gtest/gtest.h"
#include <memory>
#include <iostream>

import qchem.Hamiltonian;
import qchem.Hamiltonian.Internal.Hamiltonian;   // tHamiltonianImp -- the term-list fold under test
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Types;

using namespace qchem;
using namespace qchem::Hamiltonian;

namespace {

//! A dynamic term that records how often the phase asked it to warm.  Everything else is the cheapest
//! legal answer -- this test is about the fold, not about physics.
class SpyDynamic : public virtual rDynamic_HT
{
public:
    mutable int refreshes=0;
    mutable int matrices=0;
    mutable int slotPreps=0;
    virtual const rsmat_t& GetMatrix(const robs_t*, const Spin&, const rChargeDensity*) const override
    { ++matrices; return itsM; }
    virtual void RefreshForDensity(const rChargeDensity*) const override {++refreshes;}
    virtual void PrepareSlots(const rbs_t*) const override {++slotPreps;}
    virtual void GetEnergy(EnergyBreakdown&, const rDM_CD*) const override {}
    virtual std::ostream& Write(std::ostream& os) const override {return os;}
private:
    rsmat_t itsM;
};

//! A dynamic term whose honest answer to BOTH phases is "nothing to do".
//!
//! ⛔ IT USED TO BE "a term that does NOT override the hook", and that class no longer exists: both hooks
//! are PURE since 2026-09-10 (user).  The defaulted no-op was hiding three terms that DID need the phase --
//! two refitting lazily inside the block loop, one a forwarding term whose children the fold never reached.
//! ⇒ A term with genuinely nothing to do now says so explicitly, which is the only way the ones that DO
//! have something to do stay visible.
class SilentDynamic : public virtual rDynamic_HT
{
public:
    virtual const rsmat_t& GetMatrix(const robs_t*, const Spin&, const rChargeDensity*) const override
    { return itsM; }
    virtual void RefreshForDensity(const rChargeDensity*) const override {}   // no memo: nothing to warm
    virtual void PrepareSlots(const rbs_t*) const override {}                 // no cache: no slots
    virtual void GetEnergy(EnergyBreakdown&, const rDM_CD*) const override {}
    virtual std::ostream& Write(std::ostream& os) const override {return os;}
private:
    rsmat_t itsM;
};

//! A static term.  It has no density and therefore nothing to warm; the phase must never reach it.
class SpyStatic : public virtual rStatic_HT
{
public:
    mutable int touched=0;
    mutable int slotPreps=0;
    virtual const rsmat_t& GetMatrix(const robs_t*, const Spin&) const override {++touched; return itsM;}
    virtual void PrepareSlots(const rbs_t*) const override {++slotPreps;}
    virtual void GetEnergy(EnergyBreakdown&, const rDM_CD*) const override {}
    virtual std::ostream& Write(std::ostream& os) const override {return os;}
private:
    rsmat_t itsM;
};

} // anon

//---------------------------------------------------------------------------------------------------
TEST(EagerRefresh, TheHamiltonianFoldsTheRefreshOverEveryDynamicTermExactlyOnce)
{
    tHamiltonianImp<double> H;
    auto* a=new SpyDynamic;  auto* b=new SpyDynamic;  auto* st=new SpyStatic;
    H.Add(a); H.Add(b); H.Add(new SilentDynamic); H.Add(st);

    // A null density is the "nothing to refresh" case and must be a quiet no-op, not a crash: the SCF
    // reaches this phase on paths where the density is not yet built.
    H.RefreshForDensity(nullptr, nullptr);
    EXPECT_EQ(a->refreshes, 0);
    EXPECT_EQ(b->refreshes, 0);

    const rChargeDensity* cd=reinterpret_cast<const rChargeDensity*>(0x1);   // never dereferenced: the
    H.RefreshForDensity(nullptr, cd);                                                 // spies ignore it
    EXPECT_EQ(a->refreshes, 1) << "the phase must reach every dynamic term";
    EXPECT_EQ(b->refreshes, 1) << "the phase must reach every dynamic term";
    EXPECT_EQ(st->touched,  0) << "a STATIC term is density-independent -- the phase must not touch it";

    // Once per call, not once per block: the whole point is that N blocks cost ONE warm.
    H.RefreshForDensity(nullptr, cd);
    EXPECT_EQ(a->refreshes, 2);
    EXPECT_EQ(a->matrices,  0) << "the refresh phase must not assemble any matrix";
}

//---------------------------------------------------------------------------------------------------
// The default is a no-op, and that is load-bearing: most terms have no density-dependent memo and must not
// be forced to say so.  A term that does not override the hook must still be addable and callable.
TEST(EagerRefresh, ATermThatDoesNotOverrideTheHookIsUndisturbed)
{
    tHamiltonianImp<double> H;
    H.Add(new SilentDynamic);
    EXPECT_NO_THROW(H.RefreshForDensity(nullptr, reinterpret_cast<const rChargeDensity*>(0x1)));
}

//---------------------------------------------------------------------------------------------------
// ⚠ THE CONTRACT THE PHASE DOES **NOT** MAKE, pinned so nobody strengthens it by accident.
// RefreshForDensity is a PRE-WARM, not a replacement for the lazy fill: every memo keeps its own
// density-serial guard, and those guards remain the correctness mechanism.  A term reached with an
// unrefreshed density still has to produce the right answer -- energy evaluation and the unit tests drive
// terms outside any prologue -- so an assert of the form "the phase must have run first" would be false.
TEST(EagerRefresh, AssemblyWithoutAPriorRefreshIsLegal)
{
    tHamiltonianImp<double> H;
    auto* a=new SpyDynamic;
    H.Add(a);
    EXPECT_NO_THROW(a->GetMatrix(nullptr, Spin::Up, nullptr));
    EXPECT_EQ(a->refreshes, 0) << "assembly must not require -- nor silently trigger -- the phase";
}

//---------------------------------------------------------------------------------------
//  THE SECOND DUTY (R1.0h, 2026-09-09): the phase also PRE-CREATES this iteration's cache slots.
//
//  Same reasoning as the file header: this is PLUMBING, so pin CALL COUNTS.  What the slot pre-creation
//  buys -- that the block loop performs no map INSERTION, hence no tree mutation, hence no race between
//  blocks whose keys differ -- is not observable from outside the term (itsCache is protected, and rightly
//  so).  What IS observable, and is what rots, is whether the fold reaches every term with the block list.
TEST(EagerRefresh, ThePhasePreparesSlotsOnEveryTermIncludingStatics)
{
    tHamiltonianImp<double> H;
    auto* d1=new SpyDynamic; auto* d2=new SpyDynamic; auto* st=new SpyStatic;
    H.Add(d1); H.Add(d2); H.Add(new SilentDynamic); H.Add(st);

    // ⚠ THE SLOT DUTY IS NOT GATED ON A DENSITY, and the density duty is: slots depend on the BLOCK LIST,
    // rho does not exist yet on the paths that reach here with a null density.
    H.RefreshForDensity(nullptr, nullptr);
    EXPECT_EQ(d1->slotPreps, 1);
    EXPECT_EQ(d2->slotPreps, 1);
    EXPECT_EQ(st->slotPreps, 1) << "a STATIC term has slots too -- its cache is never cleared, so it is "
                                   "iteration ONE that would otherwise insert from inside the block loop";
    EXPECT_EQ(d1->refreshes, 0) << "no density => nothing to warm";
    EXPECT_EQ(st->touched,   0) << "and a static term is still never asked to refresh FOR A DENSITY";

    const rChargeDensity* cd=reinterpret_cast<const rChargeDensity*>(0x1);   // never dereferenced
    H.RefreshForDensity(nullptr, cd);
    EXPECT_EQ(d1->slotPreps, 2) << "the slot duty runs EVERY iteration: last iteration's blocks are stale";
    EXPECT_EQ(d1->refreshes, 1);

    // A term that overrides neither hook must be undisturbed by both -- the SilentDynamic above would have
    // thrown or crashed if the fold assumed an override.
    SUCCEED();
}
