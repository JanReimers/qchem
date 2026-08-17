// File: ElectronConfigurations/OccupationPolicy.C  The occupation policy -- WHO decides how orbitals fill.
//
// V1.11 increment 3 (design: doc/SCFStrategyPlan.md §5b, user rulings 2026-08-17).  The occupation seam's
// state and configuration, moved OFF the wave function: the WF face used to carry SetMOM / SetSmearing /
// AdoptMOMReference / ReleaseMOMReference / GetEntropyTerm as post-ctor setters and side-channels, with the
// MOM reference matrices living inside every tIrrepWF.  Now the SCFIterator owns ONE policy object (its
// occupation SLOT, built like the density mixer from SCFParams) and passes it to every fill; the WF carries
// no occupation state at all.  This also makes the seed fill's occupation EXPLICIT: Init receives the
// iterator's policy in its DEFAULT state (prescribed integer fill, kT=0, no MOM) -- which is what a seed
// fill SHOULD do, and what the old two-phase configuration did implicitly while losing charge on metals
// when anyone assumed otherwise (the D11 lesson, see tCompositeWF::FillOrbitals).
//
// HOME: qcElectronConfiguration (user ruling D6) -- beside ElectronConfiguration, whose reservoir
// PARTITION this policy's fill rules act over.  qcOrbitals sits ABOVE this library in the DAG, so the
// reference capture and overlap scores read orbitals through the OrbitalView DIP face (user ruling:
// invert, never re-point the DAG) -- TOrbitals implements it from above.
//
// SHAPE (increments 4+5 landed): the policy DECIDES -- DecideBlockFill produces the BlockFill spec (the
// §5b two-axis product: occupancy {Integer, Fermi, FermiAtMu} × ranking {priority, eShift}) and the
// orbitals execute it; HeldOccupationPolicy below is the direct minimiser's sibling (the ex-holdBlock
// bool).  Still at call sites by design: the shared-μ metal spec (the reservoir-driver migration).
module;
#include <map>
export module qchem.ElectronConfiguration.OccupationPolicy;
export import qchem.ElectronConfiguration.OrbitalView;  // the DIP view of an orbital block (see its header)
export import qchem.Symmetry.Irrep;  // Irrep (the per-block key)
import qchem.Types;              // mat_t<T>, rvec_t

export namespace qchem {

//! \brief The occupation policy for one SCF run: the run's occupation CONFIGURATION (smearing kT, MOM and
//! its penalty -- ex SetMOM/SetSmearing) plus the policy-owned per-block STATE (the MOM reference occupied
//! subspaces, the delayed-capture fill counts) and the per-fill entropy aggregate.  Owned by the
//! SCFIterator (its occupation slot); handed to every WF fill by reference.  A default-constructed policy
//! IS the seed policy: prescribed integer fill, kT=0, no MOM.
template <class T> class OccupationPolicy
{
public:
    OccupationPolicy() = default;
    virtual ~OccupationPolicy() = default;

    //! Configure for the run (the SCFIterator calls this at the top of Iterate from SCFParams -- the same
    //! slot pattern as the density mixer).  References adopted BEFORE Iterate (grid continuation) survive:
    //! configuration changes, state persists -- which is also what a staged/annealed run needs between
    //! its Iterate calls.
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

    // ---- remaining run-config queries (the reservoir driver + trace/debug still consult these; they
    //      shrink further when the driver migrates here) ----
    virtual bool   UseMOM()          const {return itsUseMOM;}
    virtual double SmearingkT()      const {return itsSmearingkT;}   //!< 0 = integer fills; >0 = Fermi (Hartree)

    // ---- MOM reference state (policy-owned; ex tIrrepWF::itsRefOccCPrime + itsFillCount) ----
    bool   HasReference(const Irrep& q) const
    { auto i=itsBlocks.find(q); return i!=itsBlocks.end() && i->second.refOccCPrime.columns()>0; }
    //! Score each of \a orbs' orbitals by the norm of its projection onto \a q's reference occupied
    //! subspace (coefficients are orthonormal-basis, metric = I); empty when no reference is held.
    rvec_t Scores(const Irrep& q, const OrbitalView<T>& orbs) const;
    //! Adopt \a from's MAJORITY-filled orbitals as \a q's fixed reference -- grid continuation
    //! (a CONVERGED foreign WF's block) or the delayed-IMOM self-capture.
    void   AdoptReference(const Irrep& q, const OrbitalView<T>& from);
    //! The delayed-IMOM capture (doc/GPWPlan §0b″): capture \a own as \a q's reference ONCE, after the
    //! settling delay -- when MOM is on, no reference is held yet, and \a q has filled >= startIter times.
    virtual void CaptureReferenceIfDue(const Irrep& q, const OrbitalView<T>& own)
    {
        auto i=itsBlocks.find(q);
        if (itsUseMOM && i!=itsBlocks.end() && i->second.refOccCPrime.columns()==0
                      && i->second.fillCount>=itsMOMStartIter) AdoptReference(q,own);
    }
    virtual void CountFill(const Irrep& q) { ++itsBlocks[q].fillCount; }   //!< one per block fill (the IMOM delay clock)
    //! The 0h MOM-guard actuator: drop EVERY reference and restart the fill counts, so the next fills run
    //! plain aufbau for the settling window and a fresh reference is captured from the recovered state.
    void   ReleaseReferences();

    // ---- cross-irrep (molecular, parked) MOM arming -- ex tCompositeWF::itsMOMActive ----
    virtual bool CrossIrrepMOMArmed() const {return itsCrossIrrepArmed;}
    void   ArmCrossIrrepMOM()         {if (itsUseMOM) itsCrossIrrepArmed=true;}

    // ---- per-fill aggregation (transitional home until increment 4's named reservoir-fill result) ----
    virtual void BeginFill()                        {itsMinusTS=0.0;}   //!< composite fill entry: reset the aggregate
    virtual void AccumulateEntropy(double wMinusTS) {itsMinusTS+=wMinusTS;}  //!< one block's BZ-weighted −TS
    //! The run's Mermin \f$-TS=\sum_k w_k(-TS_k)\le0\f$ from the most recent fill; 0 unless smearing is on.
    //! The SCFIterator stamps it into EnergyBreakdown::MinusTS (and gates direct-min on A=E−TS with it).
    virtual double EntropyTerm() const              {return itsMinusTS;}

protected:
    struct BlockState
    {
        mat_t<T> refOccCPrime;  //!< MOM reference: occupied C' columns (nbasis x nocc); empty = none
        int      fillCount=0;   //!< # fills of this block (≈ SCF iteration) -- the IMOM capture delay
    };
    std::map<Irrep,BlockState> itsBlocks;
    bool   itsUseMOM=false;
    int    itsMOMStartIter=10;
    double itsSmearingkT=0.0;
    double itsMOMSmearPenalty=0.0;
    bool   itsCrossIrrepArmed=false;  //!< the parked molecular cross-irrep MOM, armed after the 1st ranked fill
    double itsMinusTS=0.0;
};

//! \brief The direct minimiser's fill policy (V1.11 increment 5 -- the §5-flagged {occupation ×
//! direct-min} cell, NAMED): KEEP THE STORED OCCUPIED BLOCK.  A geodesic rotates a FIXED occupied block
//! and returns it as the leading columns, so only a stored-order fill reproduces the block its search
//! direction was built for; a Fermi μ solved over the whole spectrum, or a MOM overlap ranking, may
//! occupy a DIFFERENT set, making E(t) discontinuous in t and the line search meaningless (measured on
//! MnO: the best of 12 backtracks still +14.5 Ha at t=2.4e-4; doc/SymmetryUpgradePlan.md §7 step 7).
//! Entropy is identically zero -- STRUCTURALLY now (SmearingkT()==0), not by documentation -- so a held
//! leg minimises E, never A=E−TS.  WRAPS the run policy: the per-block fill COUNTS still tick on the
//! shared clocks (the delayed-IMOM capture delay keeps advancing through a direct-min leg, as before) and
//! the −TS aggregate is the run policy's (the iterator keeps reading ONE object), but no reference is
//! ever captured from a trial state and no cross-block re-decision can occur (\c HoldsStoredBlocks
//! degrades the reservoir loop to per-block fills).  The \c holdBlock bool that threaded through
//! FillOrbitals/MoveOrbitals died here.  Ready for OT+smearing: a future coupled leg swaps in a policy
//! that holds the block but keeps kT -- a new sibling, not a new bool.
template <class T> class HeldOccupationPolicy : public OccupationPolicy<T>
{
public:
    explicit HeldOccupationPolicy(OccupationPolicy<T>& run) : itsRun(run) {}
    virtual BlockFill DecideBlockFill(const Irrep&, const OrbitalView<T>&, double ne) const override
    { BlockFill f; f.budget=ne; return f; }   // stored-order integer fill == exactly the geodesic's block
    virtual bool   HoldsStoredBlocks()  const override {return true;}
    virtual bool   UseMOM()             const override {return false;}  // never rank a trial by overlap
    virtual double SmearingkT()         const override {return 0.0;}    // never re-smear a trial; −TS==0
    virtual bool   CrossIrrepMOMArmed() const override {return itsRun.CrossIrrepMOMArmed();}
    virtual void   CountFill(const Irrep& q)  override {itsRun.CountFill(q);}   // the IMOM clock still ticks
    virtual void   CaptureReferenceIfDue(const Irrep&, const OrbitalView<T>&) override {}  // never from a trial
    virtual void   BeginFill()                override {itsRun.BeginFill();}
    virtual void   AccumulateEntropy(double w) override {itsRun.AccumulateEntropy(w);}
    virtual double EntropyTerm() const        override {return itsRun.EntropyTerm();}
private:
    OccupationPolicy<T>& itsRun;   //!< the run's policy: shared clocks + the −TS aggregate live THERE
};

} // namespace qchem
