// File: ElectronConfigurations/OccupationPolicy.C  The occupation seam -- WHO decides how orbitals fill,
// and WHERE the decision's memory lives.
//
// V1.11 increment 3 (design: doc/SCFStrategyPlan.md §5b, user rulings 2026-08-17) moved the occupation
// state and configuration OFF the wave function: the WF face used to carry SetMOM / SetSmearing /
// AdoptMOMReference / ReleaseMOMReference / GetEntropyTerm as post-ctor setters and side-channels, with the
// MOM reference matrices living inside every tIrrepWF.  Now the SCFIterator owns the seam (its occupation
// SLOT, built like the density mixer from SCFParams) and passes it to every fill; the WF carries no
// occupation state at all.  This also makes the seed fill's occupation EXPLICIT: Init receives the
// iterator's policy in its DEFAULT state (prescribed integer fill, kT=0, no MOM) -- which is what a seed
// fill SHOULD do, and what the old two-phase configuration did implicitly while losing charge on metals
// when anyone assumed otherwise (the D11 lesson, see tCompositeWF::FillOrbitals).
//
// R2.21 (2026-08-17) SPLIT THE SEAM IN TWO, which is what this file now holds:
//
//   OccupationState   -- the run's persistent MEMORY: per-block MOM references + fill counts, the
//                        cross-irrep arming, the −TS aggregate.  NOT templated, and that is the point:
//                        it survives reconfiguration (grid continuation adopts references before Iterate;
//                        an annealed run calls Iterate per stage with different kT), and it is SHARED by
//                        every policy acting on the run -- including a <double> policy over a real TRIM
//                        block inside a <dcmplx> run, which is what makes MOM work on a mixed mesh.
//   OccupationPolicy  -- the per-fill DECISION, and nothing else.  Templated on the block's scalar,
//                        because a decision reads that block's orbitals.
//
// WHY THE REFERENCES ARE TYPED PER BLOCK, not per run (the defect this split repairs).  A MOM reference is
// a matrix of the BLOCK's scalar (mat_t<double> for a real TRIM block, mat_t<dcmplx> for a general-k one),
// so a run-typed home cannot store one for a block of the other scalar.  Before this split a real block in
// a complex run had nowhere to put its reference and its capture THREW ("run with forceComplex") -- turning
// "the user switched MOM on" into a mid-run failure after the basis, the seed and the first Fock were
// already paid for.  The state keys each block's reference by the BLOCK's own scalar, so both live in one
// ledger and neither is a special case.
//
// HOME: qcElectronConfiguration (user ruling D6) -- beside ElectronConfiguration, whose reservoir
// PARTITION this policy's fill rules act over.  qcOrbitals sits ABOVE this library in the DAG, so the
// reference capture and overlap scores read orbitals through the OrbitalView DIP face (user ruling:
// invert, never re-point the DAG) -- TOrbitals implements it from above.
//
// THE POLICY IS A TWO-AXIS ASSEMBLY (R2.21 shape half, 2026-08-17).  It used to be ONE concrete class
// whose Configure(useMOM, momStartIter, kT, momPenalty) flags selected among four behaviours, re-branched
// on every fill.  Now the behaviour is COMPOSED, once, at assembly:
//
//     occupancy {Integer, Fermi(kT)}   x   ranking {Bare, MOM(startIter, Λ)}
//
// The occupancy axis says HOW the budget is distributed (a count-down, or a Fermi solve at this block's own
// μ); the ranking axis says in WHAT ORDER -- the null object Bare (stored/energy order = plain aufbau) or
// MOM (overlap onto the reference occupied subspace).  The axes are genuinely independent, which is the
// point: `kT>0` is answered ONCE by which occupancy object exists, never re-asked per fill, and adding a
// third occupancy or a third ranking is a new class rather than a fifth flag combination.
//
// CONFIGURATION IS A VALUE (OccupationConfig), not a set of setters.  Two things need it: the Factory that
// assembles the concretes at the top of each Iterate, and the mixed-mesh cross arm, which rebuilds the SAME
// policy at ANOTHER SCALAR over the SAME state (a real TRIM block inside a complex run -- see the state
// note above).  Both are one call, and neither can drift from the other.
//
// (The Factory takes an OccupationConfig rather than SCFParams as the item sketched: SCFParams lives in
// qcSCFIterator, which is ABOVE this library in the link DAG.  The iterator converts -- the same direction
// as every other DIP face here.)
module;
#include <map>
#include <memory>
#include <variant>
export module qchem.ElectronConfiguration.OccupationPolicy;
export import qchem.ElectronConfiguration.OrbitalView;  // the DIP view of an orbital block (see its header)
export import qchem.Symmetry.Irrep;  // Irrep (the per-block key)
import qchem.Types;              // mat_t<T>, rvec_t

export namespace qchem {

//! \brief The occupation seam's PERSISTENT state: what the run remembers between fills, between
//! reconfigurations, and across the whole SCF.  Owned by the SCFIterator beside its policy slot.
//!
//! Deliberately NOT templated.  Two independent reasons, both load-bearing:
//!   * it must survive RECONFIGURATION (a staged/annealed run re-assembles its policy per Iterate against
//!     the SAME memory -- the constraint that originally forced the policy's Configure flags), and
//!   * one run's blocks may carry DIFFERENT scalars (a real TRIM block beside complex general-k ones),
//!     and they share one entropy ledger, one fill clock and one arming flag while each keeping a
//!     reference of its OWN scalar.
class OccupationState
{
public:
    //! Does \a q hold a usable MOM reference (present AND non-empty)?  Scalar-agnostic: the caller asking
    //! "is there a reference" does not need to know the block's type.
    bool HasReference(const Irrep& q) const;
    //! \a q's MOM reference as the BLOCK's scalar -- null when none is held.  THROWS when a reference is
    //! held under the other scalar: a block does not change type mid-run, so that is a defect, not a miss.
    template <class T> const mat_t<T>* Reference(const Irrep& q) const;
    //! Score each of \a orbs' orbitals by the norm of its projection onto \a q's reference occupied
    //! subspace (coefficients are orthonormal-basis, metric = I); empty when no reference is held.
    template <class T> rvec_t Scores(const Irrep& q, const OrbitalView<T>& orbs) const;
    //! Adopt \a from's MAJORITY-filled orbitals as \a q's fixed reference -- grid continuation
    //! (a CONVERGED foreign WF's block) or the delayed-IMOM self-capture.  Stored as \a from's scalar.
    template <class T> void AdoptReference(const Irrep& q, const OrbitalView<T>& from);
    //! The 0h MOM-guard actuator: drop EVERY reference and restart the fill counts, so the next fills run
    //! plain aufbau for the settling window and a fresh reference is captured from the recovered state.
    void ReleaseReferences();

    // ---- the fill clock (the delayed-IMOM capture delay) ----
    void CountFill(const Irrep& q) {++itsBlocks[q].fillCount;}   //!< one per block fill
    int  FillCount(const Irrep& q) const;

    // ---- cross-irrep (molecular, parked) MOM arming -- ex tCompositeWF::itsMOMActive ----
    bool CrossIrrepMOMArmed() const {return itsCrossIrrepArmed;}
    void ArmCrossIrrepMOM()         {itsCrossIrrepArmed=true;}

    // ---- per-fill aggregation ----
    void   BeginFill()                        {itsMinusTS=0.0;}   //!< composite fill entry: reset the aggregate
    void   AccumulateEntropy(double wMinusTS) {itsMinusTS+=wMinusTS;}  //!< one block's BZ-weighted −TS
    //! The run's Mermin \f$-TS=\sum_k w_k(-TS_k)\le0\f$ from the most recent fill; 0 unless smearing is on.
    //! The SCFIterator stamps it into EnergyBreakdown::MinusTS (and gates direct-min on A=E−TS with it).
    double EntropyTerm() const                {return itsMinusTS;}

private:
    //! One block's memory.  The reference is keyed by the BLOCK's scalar (see the file header): monostate
    //! = none yet.  A run with a real TRIM block and complex general-k blocks stores both alternatives in
    //! one map, which is exactly the mixed-mesh case R2.21 exists to serve.
    struct BlockState
    {
        std::variant<std::monostate, mat_t<double>, mat_t<dcmplx>> ref;
        int fillCount=0;   //!< # fills of this block (≈ SCF iteration) -- the IMOM capture delay
    };
    std::map<Irrep,BlockState> itsBlocks;
    bool   itsCrossIrrepArmed=false;  //!< the parked molecular cross-irrep MOM, armed after the 1st ranked fill
    double itsMinusTS=0.0;
};

//! \brief The occupation seam's CONFIGURATION, as a value.  What the run asked for -- never how it is
//! implemented.  The Factory turns it into objects; the mixed-mesh cross arm re-runs the Factory at another
//! scalar with the SAME value, which is why this is a struct and not a set of setters (R2.21).
struct OccupationConfig
{
    bool   useMOM       = false;  //!< rank fills by overlap onto a captured reference
    int    momStartIter = 10;     //!< the delayed-IMOM settling window, in fills of the block
    double kT           = 0.0;    //!< Fermi smearing temperature, Hartree; 0 = integer fills
    double momPenalty   = 0.0;    //!< Λ: the MOM-masked-Fermi effective-energy shift (smeared runs only)
};

//! \name The OCCUPANCY axis -- HOW a block's electron budget is distributed.
//! Scalar-FREE: distributing a budget over levels never touches an orbital coefficient.
//!@{
class OccupancyRule
{
public:
    virtual ~OccupancyRule() = default;
    virtual void   Apply(BlockFill&) const = 0;   //!< stamp the rule (and its parameters) onto the spec
    virtual double kT() const {return 0.0;}       //!< 0 for every integer rule
};

//! Prescribed integer count-down -- the aufbau/seed occupancy, and \c BlockFill's own default, so this is
//! the axis's null object as well as a member of it.
class IntegerOccupancy final : public OccupancyRule
{
public:
    virtual void Apply(BlockFill&) const override {}   // BlockFill::Rule::Integer is the default
};

//! FERMI SMEARING (doc/GPWPlan1.md 4b): solve μ per block by bisection on Σg_i f_i=ne and fill
//! fractionally.  It cures the near-gapless occupation FLAPPING (NaF Ecut=160) that MOM (an
//! occupied-subspace pin) cannot -- the frontier there IS near-degenerate, so no integer configuration is
//! stable; the fractional fill is.
class FermiOccupancy final : public OccupancyRule
{
public:
    explicit FermiOccupancy(double kT) : itsKT(kT) {}
    virtual void   Apply(BlockFill& f) const override {f.rule=BlockFill::Rule::Fermi; f.kT=itsKT;}
    virtual double kT() const override {return itsKT;}
private:
    double itsKT;
};
//!@}

//! \name The RANKING axis -- in WHAT ORDER the occupancy sees the levels.
//! Scalar-TYPED: a ranking reads this block's orbitals.
//!@{
template <class T> class RankingRule
{
public:
    virtual ~RankingRule() = default;
    //! Stamp the effective-energy ranking onto \a f, whose occupancy rule is ALREADY set (the two axes
    //! meet here: \c Integer takes a priority order, \c Fermi an ε-shift -- see MOMRanking).
    virtual void Apply(BlockFill&, const Irrep&, const OrbitalView<T>&, const OccupationState&) const {}
    //! Does this ranking consult captured references?  (The run-config query the reservoir driver asks.)
    virtual bool UsesReferences() const {return false;}
    //! The end-of-fill hook: capture a reference if this ranking wants one and it is due.
    virtual void CaptureIfDue(const Irrep&, const OrbitalView<T>&, OccupationState&) const {}
};

//! The NULL OBJECT of the ranking axis: stored/energy order, i.e. plain aufbau.  Not a "MOM off" flag --
//! a member of the axis, which is what lets the policy hold a ranking unconditionally.
template <class T> class BareRanking final : public RankingRule<T> {};

//! MOM (maximum-overlap method): rank by projection onto the reference occupied subspace.
//!
//! HOW THE TWO PATHS DIFFER, and they are DIFFERENT INSTRUMENTS rather than two strengths of one:
//!  * COLD (integer occupancy): the scores ARE the priority order, so any separation, however small,
//!    decides the fill.
//!  * MASKED FERMI (smeared occupancy, Λ>0): push low-overlap ghosts UP in effective energy
//!    (ε_i + Λ(1−s_i)²) so they stay empty BY CHARACTER, while retained high-overlap physical states smear
//!    by their TRUE energy.  Only differences worth more than Λ in ENERGY decide anything.
//!
//! HOW TO SCALE Λ (measured on MnO 2026-08-08, doc/SymmetryUpgradePlan.md §7 step 7 -- and NOT what the
//! idealised s≈0.9-vs-s≈0.1 picture suggests).  Real scores do NOT split into two clean camps: MnO's first
//! fill runs 0.95 down to 0.69 with a CUT GAP of 0.0147 between the last occupied and the first virtual.
//! So the shift separating them is Λ·[(1−s_lo)²−(1−s_hi)²] ≈ 0.0086·Λ -- 2.6 mHa at Λ=0.3, against a raw-ε
//! spread of 0.2-1 Ha.  **Λ scaled to the TIE does nothing at all** (Λ=0.3 measured indistinguishable from
//! MOM off).  What a working Λ buys is the shift on genuinely FOREIGN states (s≈0.5 ⇒ Λ(0.51)² ≈ 0.39 Ha at
//! Λ=1.5), which is what keeps a diving foreign state out.  So: scale Λ to the PHYSICAL-vs-FOREIGN score
//! gap, not to the frontier tie.  QCHEM_MOM_SCORES prints the distribution to calibrate against.
//!
//! IMOM (Initial MOM): the reference is captured ONCE, after a settling delay, from the (now physical)
//! occupied subspace, then held FIXED.  Capturing too early (the raw seed, mid-transient) locks onto
//! garbage (measured: +5 Ha states occupied, −112 Ha empty); re-capturing every iteration DRIFTS and a
//! spike corrupts the reference.
template <class T> class MOMRanking final : public RankingRule<T>
{
public:
    MOMRanking(int startIter, double penalty) : itsStartIter(startIter), itsPenalty(penalty) {}
    virtual void Apply(BlockFill&, const Irrep&, const OrbitalView<T>&, const OccupationState&) const override;
    virtual bool UsesReferences() const override {return true;}
    virtual void CaptureIfDue(const Irrep& q, const OrbitalView<T>& own, OccupationState& st) const override
    { if (!st.HasReference(q) && st.FillCount(q)>=itsStartIter) st.AdoptReference(q,own); }
private:
    int    itsStartIter;
    double itsPenalty;
};
//!@}

//! \brief The occupation POLICY for one SCF run: the per-fill DECISION, and nothing else.  All memory
//! lives in the OccupationState it is built over (R2.21); all configuration arrived as a value.
//!
//! Owned by the SCFIterator (its occupation slot); handed to every WF fill by reference.  A policy
//! assembled from a default OccupationConfig IS the seed policy: prescribed integer fill, kT=0, no MOM.
//!
//! MIXED RUNS: a policy is templated on the BLOCK's scalar, not the run's.  Every policy over one run
//! shares that run's single OccupationState, so a \c OccupationPolicy<double> serving a real TRIM block
//! inside a complex run keeps the same entropy ledger, the same fill clock and the same arming flag as its
//! complex siblings -- and captures its reference into the same map, under its own scalar.
template <class T> class OccupationPolicy
{
public:
    OccupationPolicy(OccupationState& state, const OccupationConfig& cfg)
        : itsState(state), itsConfig(cfg) {}
    virtual ~OccupationPolicy() = default;

    //! \brief THE per-block fill decision (V1.11 inc 4; the §5b two-axis product): the occupancy rule
    //! stamps how the budget is distributed, the ranking stamps in what order.  \a ne is the block's
    //! electron budget -- the reservoir distribution's output.  The orbitals EXECUTE the returned spec
    //! (\c TOrbitals::Fill); the WF decides nothing.  (Only the shared-μ metal fill still builds its spec
    //! at the call site -- the reservoir-driver migration.)
    virtual BlockFill DecideBlockFill(const Irrep& q, const OrbitalView<T>& orbs, double ne) const = 0;

    //! Does this policy pin each block's STORED occupied set -- i.e. must the reservoir loop degrade to
    //! per-block fills, with no cross-block re-decision (ranking, shared μ)?  The direct minimiser's
    //! HeldOccupationPolicy answers true; the run policy false.  (V1.11 inc 5 -- the ex-holdBlock bool.)
    virtual bool   HoldsStoredBlocks() const = 0;
    virtual bool   UseMOM()            const = 0;
    virtual double SmearingkT()        const = 0;   //!< 0 = integer fills; >0 = Fermi (Hartree)

    //! \brief The END-OF-FILL hook: this block was just filled (V1.11's two calls, collapsed).  Advances
    //! the block's clock and lets the ranking capture a reference if it wants one and it is due.
    virtual void OnBlockFilled(const Irrep& q, const OrbitalView<T>& own) = 0;

    // ---- the shared ledger: forwarded so the WF face (FillOrbitals(pol,...)) never sees the state ----
    bool   HasReference(const Irrep& q) const                    {return itsState.HasReference(q);}
    rvec_t Scores(const Irrep& q, const OrbitalView<T>& o) const  {return itsState.Scores(q,o);}
    //! The shared-μ metal path's clock tick.  Deliberately NOT OnBlockFilled: that fill is a plain energy
    //! Fermi at a μ solved on bare ε (no MOM eShift entered it), so capturing a reference FROM it would
    //! snapshot a subspace the ranking never shaped.  The clock still advances.
    void   CountFill(const Irrep& q)                             {itsState.CountFill(q);}
    bool   CrossIrrepMOMArmed() const                            {return itsState.CrossIrrepMOMArmed();}
    void   ArmCrossIrrepMOM()                                    {if (UseMOM()) itsState.ArmCrossIrrepMOM();}
    void   BeginFill()                                           {itsState.BeginFill();}
    void   AccumulateEntropy(double wMinusTS)                    {itsState.AccumulateEntropy(wMinusTS);}
    double EntropyTerm() const                                   {return itsState.EntropyTerm();}

    //! The run's memory and the run's request -- so a composite can rebuild THIS policy at ANOTHER scalar
    //! (the mixed-mesh cross arm) with one Factory call and no side channel.
    OccupationState&        State()  const {return itsState;}
    const OccupationConfig& Config() const {return itsConfig;}

protected:
    OccupationState& itsState;    //!< the run's persistent memory (not owned; outlives every policy)
    OccupationConfig itsConfig;   //!< what the run asked for (the Factory's input, kept for the cross arm)
};

//! \brief The assembled run policy: one occupancy object times one ranking object, both chosen ONCE by the
//! Factory.  This class holds no flags and takes no branches -- the four behaviours Configure used to
//! select between are four pairs of objects.
template <class T> class RunOccupationPolicy final : public OccupationPolicy<T>
{
public:
    RunOccupationPolicy(OccupationState& state, const OccupationConfig& cfg,
                        std::unique_ptr<OccupancyRule> occ, std::unique_ptr<RankingRule<T>> rank)
        : OccupationPolicy<T>(state,cfg), itsOccupancy(std::move(occ)), itsRanking(std::move(rank)) {}

    virtual BlockFill DecideBlockFill(const Irrep& q, const OrbitalView<T>& orbs, double ne) const override
    {
        BlockFill f;
        f.budget=ne;
        itsOccupancy->Apply(f);                                    // WHICH rule distributes the budget
        itsRanking->Apply(f,q,orbs,this->itsState);                // in WHAT order it sees the levels
        return f;
    }
    virtual bool   HoldsStoredBlocks() const override {return false;}
    virtual bool   UseMOM()            const override {return itsRanking->UsesReferences();}
    virtual double SmearingkT()        const override {return itsOccupancy->kT();}
    virtual void   OnBlockFilled(const Irrep& q, const OrbitalView<T>& own) override
    {
        this->itsState.CountFill(q);
        itsRanking->CaptureIfDue(q,own,this->itsState);
    }
private:
    std::unique_ptr<OccupancyRule>   itsOccupancy;
    std::unique_ptr<RankingRule<T>>  itsRanking;
};

//! \brief The direct minimiser's fill policy (V1.11 increment 5 -- the §5-flagged {occupation ×
//! direct-min} cell, NAMED): KEEP THE STORED OCCUPIED BLOCK.  A geodesic rotates a FIXED occupied block
//! and returns it as the leading columns, so only a stored-order fill reproduces the block its search
//! direction was built for; a Fermi μ solved over the whole spectrum, or a MOM overlap ranking, may
//! occupy a DIFFERENT set, making E(t) discontinuous in t and the line search meaningless (measured on
//! MnO: the best of 12 backtracks still +14.5 Ha at t=2.4e-4; doc/SymmetryUpgradePlan.md §7 step 7).
//! Entropy is identically zero -- STRUCTURALLY (SmearingkT()==0), not by documentation -- so a held leg
//! minimises E, never A=E−TS.
//!
//! It is Integer×Bare with one thing added -- \c HoldsStoredBlocks -- so it is a genuine sibling of the
//! assembled policy rather than a wrapper: the shared clocks and the −TS aggregate are the STATE's, so a
//! held policy built over the same state gets them for free (the delayed-IMOM delay keeps advancing
//! through a direct-min leg, as before, and the iterator still reads ONE −TS), and no reference is ever
//! captured from a trial state.  Ready for OT+smearing: a coupled leg that holds the block but keeps kT is
//! a new sibling with a Fermi occupancy -- a new object, not a new bool.
template <class T> class HeldOccupationPolicy final : public OccupationPolicy<T>
{
public:
    explicit HeldOccupationPolicy(OccupationPolicy<T>& run) : OccupationPolicy<T>(run.State(), OccupationConfig{}) {}
    virtual BlockFill DecideBlockFill(const Irrep&, const OrbitalView<T>&, double ne) const override
    { BlockFill f; f.budget=ne; return f; }   // stored-order integer fill == exactly the geodesic's block
    virtual bool   HoldsStoredBlocks() const override {return true;}
    virtual bool   UseMOM()            const override {return false;}  // never rank a trial by overlap
    virtual double SmearingkT()        const override {return 0.0;}    // never re-smear a trial; −TS==0
    //! The clock still ticks; the capture never fires (never from a trial state).
    virtual void   OnBlockFilled(const Irrep& q, const OrbitalView<T>&) override {this->itsState.CountFill(q);}
};

//! \brief Assemble the run's policy from what it asked for.  THE one place the two axes are chosen, and it
//! runs once per Iterate -- so `kT>0` is a question about which object exists, never a per-fill branch.
template <class T> std::unique_ptr<OccupationPolicy<T>> MakeOccupationPolicy(const OccupationConfig&,
                                                                             OccupationState&);

} // namespace qchem
