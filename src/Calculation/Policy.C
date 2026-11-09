// File: Calculation/Policy.C  How a policy reports itself: one row in the run report where it fired (doc/OpenWork.md, "Every existing
// policy reports itself under policy").
//
// Policies live only in the policy layer (the facades, their setups and factories; CLAUDE.md framework-vs-efficiency), and each one
// that ACTS writes one row, at its report cursor -- the caller set the context, so the row sits where the policy fired (the run's own
// level, basis, ...), and a reader collects every `policy` array in the run's tree (ruled 2026-11-07).  The row:
//     {name, decided, verdict: refuse | resource | warn | choose, reason, cost?: {wallSeconds?, memoryBytes?}, detail?: {...}}
// `name` is the policy's identifier in the code, so grepping the report's names finds the code.  `choose` is an AUTOMATIC setting a
// policy picked because the caller stated none.  A refusal writes its row and then throws: the driver has opened the run before
// any refusal can fire, so the record of a refused run carries it.  Outside a run the row is inert (the report's own rule).
module;
#include <optional>
#include <stdexcept>
#include <string>
export module qchem.Policy;
import qchem.Reporting;

export namespace qchem::policy
{

enum class Verdict { Refuse, Resource, Warn, Choose };

//! What a warn-and-run expects to spend, or what a resource refusal would have needed; either half alone is fine.
struct Cost
{
    std::optional<double> wallSeconds;
    std::optional<double> memoryBytes;
};

inline const char* NameOf(Verdict v)
{
    switch (v)
    {
    case Verdict::Refuse:   return "refuse";
    case Verdict::Resource: return "resource";
    case Verdict::Warn:     return "warn";
    case Verdict::Choose:   return "choose";
    }
    return "?";
}

//! One decision, as a row of the `policy` array at the report cursor.  \a decided is one line (what), \a reason why.  \a detail (optional)
//! carries the numbers the decision was made from, as the policy's own json (a resource policy: its estimate, calibration and the memory
//! available), so a reader need not parse them out of \a reason.
inline void Report(const std::string& name, const std::string& decided, Verdict verdict, const std::string& reason,
                   const std::optional<Cost>& cost = std::nullopt, const qchem::report::json& detail = qchem::report::json())
{
    namespace rpt = qchem::report;
    const rpt::Row row("policy");
    rpt::Set("name",    name);
    rpt::Set("decided", decided);
    rpt::Set("verdict", NameOf(verdict));
    rpt::Set("reason",  reason);
    if (cost && (cost->wallSeconds || cost->memoryBytes))
    {
        rpt::json c = rpt::json::object();
        if (cost->wallSeconds) c["wallSeconds"] = *cost->wallSeconds;
        if (cost->memoryBytes) c["memoryBytes"] = *cost->memoryBytes;
        rpt::Set("cost", c);
    }
    if (!detail.is_null()) rpt::Set("detail", detail);
}

//! A refusal: its row, then std::invalid_argument("<decided> -- <reason>") -- a user-reachable choice the caller cannot run as asked
//! (an Outcome at a facade's Converge door; a throw from a setup, which is how the setups refuse today).
[[noreturn]] inline void Refuse(const std::string& name, const std::string& decided, const std::string& reason)
{
    Report(name, decided, Verdict::Refuse, reason);
    throw std::invalid_argument(decided + " -- " + reason);
}

} // namespace qchem::policy
