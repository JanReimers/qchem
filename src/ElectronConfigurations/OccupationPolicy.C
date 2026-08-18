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
// REMAINING HALF of R2.21 (doc/CleanupCandidates.md): the POLICY is not yet the abstract interface the
// design ruled -- Configure's flags still select the behaviour per fill (four policies in one object,
// branched inside DecideBlockFill).  The finishing move is the two-axis assembly: occupancy
// {Integer, Fermi(kT)} composed with a ranking {Bare (the null object), MOM(startIter, Λ)}, built by a
// Factory(SCFParams, state&) at the top of each Iterate so the kT branch runs ONCE at assembly instead of
// per fill; Configure then dies.  That is SHAPE work -- the STATE half above is what the mixed-mesh MOM
// case needed, and it is done.  When it lands, note the cross-scalar rebuild this file already relies on
// (CopyConfigTo) becomes Factory<double>(config, state) -- so keep the configuration a value, not a set
// of flags smeared over the concretes.
module;
#include <map>
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

//! \brief The occupation POLICY for one SCF run: the run's occupation CONFIGURATION (smearing kT, MOM and
//! its penalty -- ex SetMOM/SetSmearing) and the per-fill DECISION it implies.  All memory lives in the
//! OccupationState this policy is built over (R2.21); the policy itself is a decision function.
//!
//! Owned by the SCFIterator (its occupation slot); handed to every WF fill by reference.  A
//! default-configured policy IS the seed policy: prescribed integer fill, kT=0, no MOM.
//!
//! MIXED RUNS: a policy is templated on the BLOCK's scalar, not the run's.  Every policy over one run
//! shares that run's single OccupationState, so a \c OccupationPolicy<double> serving a real TRIM block
//! inside a complex run keeps the same entropy ledger, the same fill clock and the same arming flag as its
//! complex siblings -- and captures its reference into the same map, under its own scalar.
template <class T> class OccupationPolicy
{
public:
    explicit OccupationPolicy(OccupationState& state) : itsState(state) {}
    virtual ~OccupationPolicy() = default;

    //! Configure for the run (the SCFIterator calls this at the top of Iterate from SCFParams -- the same
    //! slot pattern as the density mixer).  References adopted BEFORE Iterate (grid continuation) survive:
    //! configuration changes, STATE persists -- structurally now, since the state is a separate object.
    void Configure(bool useMOM, int momStartIter, double kT, double momPenalty)
    { itsUseMOM=useMOM; itsMOMStartIter=momStartIter; itsSmearingkT=kT; itsMOMSmearPenalty=momPenalty; }

    //! \brief THE per-block fill decision (V1.11 inc 4; the §5b two-axis product): the occupancy rule from
    //! the run configuration (kT), the effective-energy ranking from this block's MOM state.  \a ne is the
    //! block's electron budget -- the reservoir distribution's output.  The orbitals EXECUTE the returned
    //! spec (\c TOrbitals::Fill); the WF no longer decides anything.  (Only the shared-μ metal fill still
    //! builds its spec at the call site -- the reservoir-driver migration.)
    virtual BlockFill DecideBlockFill(const Irrep& q, const OrbitalView<T>& orbs, double ne) const;

    //! Does this policy pin each block's STORED occupied set -- i.e. must the reservoir loop degrade to
    //! per-block fills, with no cross-block re-decision (ranking, shared μ)?  The direct minimiser's
    //! HeldOccupationPolicy answers true; the run policy false.  (V1.11 inc 5 -- the ex-holdBlock bool.)
    virtual bool HoldsStoredBlocks() const {return false;}

    // ---- run-config queries (the reservoir driver + trace/debug still consult these) ----
    virtual bool   UseMOM()          const {return itsUseMOM;}
    virtual double SmearingkT()      const {return itsSmearingkT;}   //!< 0 = integer fills; >0 = Fermi (Hartree)

    //! The delayed-IMOM capture (doc/GPWPlan §0b″): capture \a own as \a q's reference ONCE, after the
    //! settling delay -- when MOM is on, no reference is held yet, and \a q has filled >= startIter times.
    virtual void CaptureReferenceIfDue(const Irrep& q, const OrbitalView<T>& own)
    {
        if (UseMOM() && !itsState.HasReference(q) && itsState.FillCount(q)>=itsMOMStartIter)
            itsState.AdoptReference(q,own);
    }

    // ---- the shared ledger: forwarded so the WF face (FillOrbitals(pol,...)) never sees the state ----
    bool   HasReference(const Irrep& q) const              {return itsState.HasReference(q);}
    rvec_t Scores(const Irrep& q, const OrbitalView<T>& o) const {return itsState.Scores(q,o);}
    virtual void   CountFill(const Irrep& q)               {itsState.CountFill(q);}
    virtual bool   CrossIrrepMOMArmed() const              {return itsState.CrossIrrepMOMArmed();}
    void           ArmCrossIrrepMOM()                      {if (itsUseMOM) itsState.ArmCrossIrrepMOM();}
    virtual void   BeginFill()                             {itsState.BeginFill();}
    virtual void   AccumulateEntropy(double wMinusTS)      {itsState.AccumulateEntropy(wMinusTS);}
    virtual double EntropyTerm() const                     {return itsState.EntropyTerm();}

    //! The run's memory -- so a composite can build a policy of ANOTHER scalar over the SAME state (the
    //! mixed-mesh cross arm) without a side channel.
    OccupationState& State() const {return itsState;}
    //! Reproduce this policy's DECISION configuration on \a p, a policy of another scalar over the same
    //! state -- the mixed-mesh cross arm's one call.  Read through the VIRTUALS, not the fields, so a
    //! policy that overrides its answers (HeldOccupationPolicy: no MOM, no smearing) is reproduced by what
    //! it ANSWERS rather than by what it stores; a held leg's real block then decides exactly as its
    //! complex siblings do.
    void CopyConfigTo(OccupationPolicy<double>& p) const
    { p.Configure(UseMOM(), itsMOMStartIter, SmearingkT(), itsMOMSmearPenalty); }

protected:
    OccupationState& itsState;        //!< the run's persistent memory (not owned; outlives every policy)
    bool   itsUseMOM=false;
    int    itsMOMStartIter=10;
    double itsSmearingkT=0.0;
    double itsMOMSmearPenalty=0.0;
};

//! \brief The direct minimiser's fill policy (V1.11 increment 5 -- the §5-flagged {occupation ×
//! direct-min} cell, NAMED): KEEP THE STORED OCCUPIED BLOCK.  A geodesic rotates a FIXED occupied block
//! and returns it as the leading columns, so only a stored-order fill reproduces the block its search
//! direction was built for; a Fermi μ solved over the whole spectrum, or a MOM overlap ranking, may
//! occupy a DIFFERENT set, making E(t) discontinuous in t and the line search meaningless (measured on
//! MnO: the best of 12 backtracks still +14.5 Ha at t=2.4e-4; doc/SymmetryUpgradePlan.md §7 step 7).
//! Entropy is identically zero -- STRUCTURALLY (SmearingkT()==0), not by documentation -- so a held
//! leg minimises E, never A=E−TS.
//!
//! R2.21: it no longer WRAPS the run policy.  The shared clocks and the −TS aggregate are the STATE's, so
//! a held policy built over the same state gets them for free -- the per-block fill counts still tick on
//! the shared clock (the delayed-IMOM delay keeps advancing through a direct-min leg, as before) and the
//! iterator still reads ONE −TS.  What remains here is exactly the three DECISIONS that differ.  The
//! \c holdBlock bool that threaded through FillOrbitals/MoveOrbitals died in V1.11.  Ready for OT+smearing:
//! a future coupled leg swaps in a policy that holds the block but keeps kT -- a new sibling, not a new bool.
template <class T> class HeldOccupationPolicy : public OccupationPolicy<T>
{
public:
    explicit HeldOccupationPolicy(OccupationPolicy<T>& run) : OccupationPolicy<T>(run.State()) {}
    virtual BlockFill DecideBlockFill(const Irrep&, const OrbitalView<T>&, double ne) const override
    { BlockFill f; f.budget=ne; return f; }   // stored-order integer fill == exactly the geodesic's block
    virtual bool   HoldsStoredBlocks()  const override {return true;}
    virtual bool   UseMOM()             const override {return false;}  // never rank a trial by overlap
    virtual double SmearingkT()         const override {return 0.0;}    // never re-smear a trial; −TS==0
    virtual void   CaptureReferenceIfDue(const Irrep&, const OrbitalView<T>&) override {}  // never from a trial
};

} // namespace qchem
