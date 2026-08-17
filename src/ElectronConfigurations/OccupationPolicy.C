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
// TRANSITIONAL SHAPE: at this increment the policy is a CONCRETE state/config carrier and the fill
// MECHANICS still live in tIrrepWF/tCompositeWF, consulting it.  Increment 4 migrates the mechanics here
// and the §5b two-axis form emerges (occupancy {Integer, Fermi, Held} × ranking {bare, MOM}), with the
// reservoir driver in the base; increment 5 replaces the holdBlock bool with the Held policy.
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
    //! spec (\c TOrbitals::Fill); the WF no longer decides anything.  (The held direct-min fill and the
    //! shared-μ metal fill still build their specs at the call site -- increment 5 / the reservoir-driver
    //! migration respectively.)
    BlockFill DecideBlockFill(const Irrep& q, const OrbitalView<T>& orbs, double ne) const;

    // ---- remaining run-config queries (the reservoir driver + trace/debug still consult these; they
    //      shrink further when the driver migrates here) ----
    bool   UseMOM()          const {return itsUseMOM;}
    double SmearingkT()      const {return itsSmearingkT;}   //!< 0 = integer fills; >0 = Fermi (Hartree)

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
    void   CaptureReferenceIfDue(const Irrep& q, const OrbitalView<T>& own)
    {
        auto i=itsBlocks.find(q);
        if (itsUseMOM && i!=itsBlocks.end() && i->second.refOccCPrime.columns()==0
                      && i->second.fillCount>=itsMOMStartIter) AdoptReference(q,own);
    }
    void   CountFill(const Irrep& q) { ++itsBlocks[q].fillCount; }   //!< one per block fill (the IMOM delay clock)
    //! The 0h MOM-guard actuator: drop EVERY reference and restart the fill counts, so the next fills run
    //! plain aufbau for the settling window and a fresh reference is captured from the recovered state.
    void   ReleaseReferences();

    // ---- cross-irrep (molecular, parked) MOM arming -- ex tCompositeWF::itsMOMActive ----
    bool   CrossIrrepMOMArmed() const {return itsCrossIrrepArmed;}
    void   ArmCrossIrrepMOM()         {if (itsUseMOM) itsCrossIrrepArmed=true;}

    // ---- per-fill aggregation (transitional home until increment 4's named reservoir-fill result) ----
    void   BeginFill()                        {itsMinusTS=0.0;}   //!< composite fill entry: reset the aggregate
    void   AccumulateEntropy(double wMinusTS) {itsMinusTS+=wMinusTS;}  //!< one block's BZ-weighted −TS
    //! The run's Mermin \f$-TS=\sum_k w_k(-TS_k)\le0\f$ from the most recent fill; 0 unless smearing is on.
    //! The SCFIterator stamps it into EnergyBreakdown::MinusTS (and gates direct-min on A=E−TS with it).
    double EntropyTerm() const                {return itsMinusTS;}

private:
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

} // namespace qchem
