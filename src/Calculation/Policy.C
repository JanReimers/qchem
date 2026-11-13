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
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
export module qchem.Policy;
import qchem.Reporting;
import qchem.Hamiltonian.Factory;   // PPFunctionalCheck (PPFittedFunctional's input, the pre-run query's record)

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

//! \brief A PSEUDOPOTENTIAL FITTED WITH ANOTHER FUNCTIONAL THAN THE RUN'S: WARN AND GUIDE, AND RUN (user rulings 2026-11-12: OpenWork PADE
//! row option a; and "most users ... will not know PPs well enough to be aware they are tuned for a particular exchange-correlation
//! functional ... the lib code should be independently user friendly and spit warnings as well").  A model choice the code can SEE
//! (pin 31: exempt from eps, so it warns and never substitutes): GTH-PADE under Dirac + VWN5 is ~0.6 mHa on the Li q3 atom, all E_xc.
//! \a checks: \c Hamiltonian::PPFunctionalChecks, the same records a GUI reads before the run; the comparison is on libxc ids, the CP2K
//! name only labels it.  Silent when every element matches; else ONE row at the cursor naming the run's functional, every mismatched
//! element with its fit (label and ids), and what to do; the same on the console.
inline void PPFittedFunctional(const std::vector<Hamiltonian::PPFunctionalCheck>& checks)
{
    std::vector<Hamiltonian::PPFunctionalCheck> mismatches;
    for (const auto& c : checks) if (!c.match) mismatches.push_back(c);
    if (mismatches.empty()) return;
    qchem::report::json detail=qchem::report::json::object();
    std::string list, advice;
    std::vector<std::string> advised;   // one piece of advice per distinct fit
    for (const auto& m : mismatches)
    {
        detail[m.element]={{"fitted",m.fittedLabel},{"fittedLibxcIDs",m.fittedIDs},{"runLibxcIDs",m.runIDs},{"matchingModel",m.matchingModel}};
        list+=(list.empty() ? "" : ", ")+m.element+" (fitted with "+m.fittedLabel+")";
        bool seen=false;
        for (const auto& f : advised) seen = seen || f==m.fittedLabel;
        if (seen) continue;
        advised.push_back(m.fittedLabel);
        advice+=(advice.empty() ? "" : "; ")+(m.matchingModel.empty()
            ? "no model offered matches "+m.fittedLabel+": use pseudopotentials fitted for "+m.runLabel
            : "state model "+m.matchingModel+" to match the GTH-"+m.fittedLabel+" fit, or use pseudopotentials fitted for "+m.runLabel);
    }
    const std::string& run=mismatches.front().runLabel;
    const std::string decided="run "+run+" with pseudopotentials fitted with another functional: "+list;
    const std::string reason ="a GTH pseudopotential is fitted with one exchange-correlation functional and is most accurate under it; this run "
                              "states "+run+" (Li q3: PADE vs Dirac + VWN5 is 0.6 mHa, all of it E_xc).  The answer is the stated model's, "
                              "solved to eps.  To match: "+advice;
    Report("PPFittedFunctional", decided, Verdict::Warn, reason, std::nullopt, detail);
    std::cout << "[policy] WARN PPFittedFunctional: " << decided << " -- " << advice << std::endl;
}

} // namespace qchem::policy
