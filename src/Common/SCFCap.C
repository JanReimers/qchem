// File: SCFCap.C  An SCF that ENDS AT ITS ITERATION CAP unconverged is reported -- and a test main makes that a failure.
//
// User ruling 2026-10-29: "All SCF tests should fail if they hit max iterations."  The A_DFT_atom suite had measured the
// energy of iteration 20 of runs that never converged, for months, and nothing said so.  A per-test assert gets forgotten,
// so the rule is enforced ONCE, at the main level: the SCF reports every cap hit here, and the two gtest mains install a
// hook that turns the report into ADD_FAILURE inside the running test (same thread -- the SCF runs on the test's own).
// No hook is installed in a production driver: there the facade's own Outcome / diagnostics speak.
//
// A test that DELIBERATELY runs a bounded SCF -- the loop is its subject, not the fixed point -- says so for its scope with
// a BoundedSCF token naming why.  Rare by design, and greppable.
//
// THE CONVERGENCE PANEL (doc/SCFRobustnessPlan.md §4, 2026-10-29) rides the same hook: a panel case holds an IterationBudget,
// the SCF reports the iterations of every Iterate it ends, and a case over its budget fails like a cap hit -- so a recipe
// change that helps one hard case and slows another goes red on the panel, not in a campaign.
module;
#include <cstddef>
#include <functional>
#include <string>
export module qchem.SCFCap;

export namespace qchem
{

//! What a cap hit is reported to: \a what names the run (its label, the cap, the last residuals).
using SCFCapHook = std::function<void(const std::string& what)>;
//! Install (or clear, with an empty function) the process-wide hook.  The test mains install ADD_FAILURE.
void SetSCFCapHook(SCFCapHook hook);
//! Called by the SCF when an Iterate ends at NMaxIter WITHOUT converging.  Suppressed inside a BoundedSCF scope.
void ReportSCFCapHit(const std::string& what);

//! Called by the SCF at the end of EVERY Iterate with the iterations it ran; every live IterationBudget on this thread counts them.
void ReportSCFIterations(std::size_t n);

//! \brief A convergence-panel budget for this scope: the SCF iterations of every Iterate (every stage) run inside it may not
//! exceed \a budget.  The destructor prints `[panel] <what>: N SCF iterations (budget B)` (the panel's harvestable line) and,
//! over budget, reports through the cap hook -- a test failure.  \a what names the case and its row or recipe.  Scopes nest
//! (an outer budget counts what the inner ones count).  A BoundedSCF does not silence it.
class IterationBudget
{
public:
    IterationBudget(std::size_t budget, const char* what);
    ~IterationBudget();
    IterationBudget(const IterationBudget&)=delete;
    IterationBudget& operator=(const IterationBudget&)=delete;
    std::size_t Used  () const {return itsUsed;}
    std::size_t Audit () const {return itsAudit;}   //!< iterations run inside an AuditSCF scope: printed, never budgeted
    std::size_t Budget() const {return itsBudget;}
private:
    friend void ReportSCFIterations(std::size_t);
    std::size_t      itsBudget;
    std::size_t      itsUsed=0;
    std::size_t      itsAudit=0;
    const char*      itsWhat;
    IterationBudget* itsOuter;   //!< the enclosing budget on this thread (a stack)
};

//! \brief A DELIBERATELY bounded SCF, for this scope: cap hits are expected and not reported.  \a why is required -- say
//! what the bounded loop is for (e.g. "the gate tests symmetry through 6 iterations, not convergence").
class BoundedSCF
{
public:
    explicit BoundedSCF(const char* why);
    ~BoundedSCF();
    BoundedSCF(const BoundedSCF&)=delete;
    BoundedSCF& operator=(const BoundedSCF&)=delete;
};

//! \brief SETUP SCFs, for this scope: the iterations of an Iterate run inside are not a panel case's to count (no live
//! IterationBudget sees them).  The facade wraps its atomic solves in one (the seed's, and the +U projector's pseudo-atom) -- a
//! budget measures the run's convergence, and the pseudo-atom solves added 122 iterations to every beta-MnO2 case (2026-10-29).
class SetupSCF
{
public:
    SetupSCF();
    ~SetupSCF();
    SetupSCF(const SetupSCF&)=delete;
    SetupSCF& operator=(const SetupSCF&)=delete;
private:
    IterationBudget* itsSaved;   //!< the budget stack, restored on exit
};

//! \brief AUDIT SCFs, for this scope: the re-convergence an audit runs AFTER the run converged (the raster pilot's, when a raster
//! rises).  A live IterationBudget counts them APART from the run's own iterations and prints them, but never gates on them: a
//! budget gates the SCF recipe, and the audit's cost is the grid's question (the pilot reports it under `rasterPilot`).  Earned by
//! the TiO2 k222 panel case (OpenWork "SCF panel iteration budgets predate the raster pilot"): its SCF converged in 43 of 49, then
//! the pilot found Hartree under-resolved and its 17 re-convergence iterations failed the case.
class AuditSCF
{
public:
    AuditSCF();
    ~AuditSCF();
    AuditSCF(const AuditSCF&)=delete;
    AuditSCF& operator=(const AuditSCF&)=delete;
};

} // namespace
