// File: Common/RunPolicy.C  THE ONE PLACE A RUN'S DEVIATIONS FROM CP2K ARE DECIDED AND NAMED.
//
// WHY THIS EXISTS (doc/OpenWork.md N5, user 2026-08-25: "there will eventually be a number of CP2K
// deviations so we will need one env bool (like CP2K_COMPAT) to trigger off").  This tree runs several
// accelerations and route choices that CP2K does not, and until now each one lived as its own
// `std::getenv` at the point of use.  Two things follow from that, and both bit:
//
//   1. NO RUN COULD STATE WHAT IT WAS.  doc/Benchmark.md's rows are only comparable if every row
//      declares its accelerations, and "declare them" was a matter of DISCIPLINE -- so the low-rank rho
//      route has been ON BY DEFAULT since 07d13bf6 and no timing row since says so.
//   2. NO RUN COULD TURN THEM ALL OFF.  Reproducing a CP2K number meant remembering six env vars.
//
// THE DESIGN, and the reason it is a policy OBJECT rather than a tidier set of getenvs: a flag read deep
// inside a kernel cannot be reported, cannot be overridden coherently, and cannot be reasoned about --
// so the choices are RESOLVED ONCE, here, and consulted where the choice is MADE (the factories),
// exactly as SCFParams::XCCuspDeficit already is.  That is what makes "CP2K parity is a property of
// what was BUILT" true rather than aspirational, and it is what lets a banner print the RESOLVED state
// instead of re-deriving it and hoping the two agree.
//
// THE OVERRIDE RULE: an explicitly-set individual knob WINS over CP2K_COMPAT.  Saying
// `CP2K_COMPAT=1 GPW_STREAM_FOLD=1` is a deliberate act -- "parity everywhere except this one" -- and
// silently overruling it would make the umbrella a trap.  The banner marks which is which, so a run
// that thinks it has parity and does not says so out loud.
//
// STANDING RULE (doc/Benchmark.md): a new accelerator is NOT FINISHED until it is in the table below.
module;
#include <cstdlib>   // std::getenv/std::atoi
#include <string>
#include <sstream>
#include <vector>
export module qchem.RunPolicy;

export namespace qchem
{

//! \brief ONE deviation: a named route or acceleration this tree runs and CP2K does not.
//! Carried as data (not as prose in a comment) so the banner is generated from the SAME facts the
//! factories consult -- there is no second list to keep in step.
struct Deviation
{
    const char* knob;    //!< the env var a user would set, e.g. "QCHEM_DM_LOWRANK"
    const char* what;    //!< one line: what it changes
    bool cp2kValue;      //!< the value CP2K parity requires
    bool value;          //!< what THIS process resolved to
    bool stated;         //!< the user named this knob explicitly (so CP2K_COMPAT did not overrule it)
    bool Deviates() const {return value!=cp2kValue;}
};

//! \brief The process-wide resolved answer to "which non-CP2K routes is this run taking?".
//!
//! Read ONCE (the values are process lifetime constants -- they select which objects the factories
//! BUILD, and re-reading them mid-run would mean two halves of one SCF disagreeing about what it is).
class RunPolicy
{
public:
    RunPolicy();

    //! \name The routes.  Each is consulted by ONE factory; see the .C file's site list.
    //!@{
    //! Factored / low-rank \f$\rho\f$ (a singles route CP2K does not run).  DEFAULT ON.
    bool DMLowRank()  const {return itsDMLowRank.value;}
    //! The T3 pair-stream orbit fold in the GPW collocation streams.  DEFAULT ON (armed 2026-08-19).
    bool StreamFold() const {return itsStreamFold.value;}
    //! Mix in the \f$(\rho,m)\f$ channel basis (Kerker on \f$\rho\f$, plain linear on \f$m\f$) instead
    //! of \f$(\rho_\uparrow,\rho_\downarrow)\f$.  DEFAULT OFF -- (up,dn) is what reproduces CP2K's Kerker.
    bool MixRhoM()    const {return itsMixRhoM.value;}
    //! Feed \f$V_{xc}\f$ the WHOLESALE \f$\rho[D]\f$ instead of Hartree's own mixed array.  DEFAULT OFF
    //! (doc/OpenWork.md item 1: measured, and the route -- not the goal -- is what fails).
    bool XCFromDM()   const {return itsXCFromDM.value;}
    //! \brief May a run EXPLOIT SPACE-GROUP SYMMETRY at all?  DEFAULT true (obey the caller); false under
    //! \c CP2K_COMPAT, which forces \c SolidCalcOptions::imposeSymmetry off however the caller set it.
    //!
    //! WHY THIS ONE OVERRULES THE CALLER, where the others merely have defaults (user, 2026-08-26:
    //! *"CP2K_COMPAT should do (imply) imposeSymmetry=0"*).  Verified from the CP2K logs: QuickStep does
    //! NO symmetry work in the decks these benchmarks compare against -- no irrep/SALC blocking (K and P
    //! are DBCSR sparse ATOM-BLOCK matrices over the full AO basis), and BZ point-group symmetrization
    //! defaults OFF.  An imposed qchem run folds the BZ, star-averages \f$\rho\f$ every iteration, uses
    //! the site-adapted invariant XC mesh and folds the collocation streams.  So leaving the imposition
    //! ON under a switch called "compat" would keep the single largest advantage the comparison is
    //! supposed to remove, and every banked recipe sets it -- a switch that needed the recipe edited too
    //! would not be one switch.  \c QCHEM_IMPOSE_SYMMETRY=1 is the stated escape hatch.
    bool SymmetryImposition() const {return itsImpose.value;}
    //! \brief May a run use the ATOM-CENTRED (Becke) XC quadrature?  DEFAULT true; false under
    //! \c CP2K_COMPAT, which then routes XC onto the uniform grid however the caller set \c xcMesh.
    //!
    //! WHY IT IS ON THIS LIST (2026-08-28).  CP2K evaluates \f$V_{xc}\f$ on the uniform realspace grid;
    //! qchem's default since 2026-08-02 is a periodic Becke atom-centred mesh, and on the MnO benchmark row
    //! that mesh is **158 s of 367 s wall — 43%, the single largest block**, ahead of the collocation.  A
    //! head-to-head row that leaves it on is qchem-with-Becke against CP2K-plain.  It reached this table
    //! late because it is a TYPED option (\c SolidCalcOptions::xcMesh) rather than an env flag — the same
    //! gap \c raster and \c cutoffFactor still sit in.
    //! ⚠ Like \c SymmetryImposition this OVERRULES the caller, for the same reason: every banked recipe
    //! sets the mesh, and a switch that needed the recipe edited too would not be one switch.
    bool BeckeXC() const {return itsBeckeXC.value;}
    //!@}

    bool CP2KCompat() const {return itsCP2KCompat;}   //!< the umbrella was asked for
    //! Every deviation, in one list, whatever its value -- the banner prints the WHOLE table, because a
    //! row that lists only what is ON cannot be read as evidence that the rest is OFF.
    std::vector<Deviation> Deviations() const
    {return {itsDMLowRank, itsStreamFold, itsMixRhoM, itsXCFromDM, itsImpose, itsBeckeXC};}
    //! Are we actually at parity?  (CP2K_COMPAT=1 plus an explicit knob that contradicts it is NOT.)
    bool AtParity() const;
    //! One line naming every deviation and its state, with `*` on the ones that differ from CP2K.
    std::string Banner() const;

private:
    Deviation Resolve(const char* knob, const char* what, bool cp2kValue, bool qchemDefault);
    bool      itsCP2KCompat = false;
    Deviation itsDMLowRank{}, itsStreamFold{}, itsMixRhoM{}, itsXCFromDM{}, itsImpose{}, itsBeckeXC{};
};

//! The process's policy.  A function rather than a global so it is constructed on first use, after
//! main() has had its chance to set the environment.
const RunPolicy& theRunPolicy();

//! \brief RE-READ the environment into the process policy.
//!
//! For A/B GATES ONLY.  Two acceptance tests flip \c GPW_STREAM_FOLD with \c setenv and run both arms
//! in ONE process, which is exactly the kind of thing a resolved-once policy forbids -- so the escape
//! hatch is named, rather than the policy being silently re-read on every call (which would make
//! "resolved once" a comment instead of a property).  A production run never calls this: half an SCF
//! disagreeing with the other half about which routes it is taking is not a state worth supporting.
void ReresolveRunPolicy();

} //export namespace qchem
