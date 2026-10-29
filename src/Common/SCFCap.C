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
module;
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

} // namespace
