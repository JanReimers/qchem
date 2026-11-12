// File: Calculation/SCFRecipeReport.C  The run report's `scf` header -- "what recipe am I running" -- for every facade.
//
// Its own module (V-FACADE-PARITY, 2026-11-02): the molecular and atomic facades each built their own header (two copies of the
// same knobs, the atom's a subset), and the solid facade filed none.  The header is a function of the SCF KNOBS, which are
// structure-neutral; what is structure-specific (the atom's model / spin / electron count) the facade hands in as its own keys.
module;
#include <iostream>
#include <sstream>
#include <string>
export module qchem.SCFRecipeReport;
export import qchem.SCFParams;
export import qchem.Reporting;
export import qchem.Accuracy;   // AccuracyTarget (ReportAccuracy's argument)
import qchem.Parallel;          // WorkerThreads, BlasThreads: the threads the recipe runs on

export namespace qchem
{

//! \brief Write the SCF recipe ONE LEVEL IN (THE RUN REPORT: the caller opens the `scf` section, this writes into it).
//!
//! \c standard = the settled recipe surface, \c advanced = the pathological-case dials (doc/RunReportPlan.md schema).  The
//! mixer is DERIVED (the preconditioner and the history COMPOSE, so both are named: "Ker+Pul", "Lin+Pul"); \a accelerator is
//! the engine tag the facade chose.  \a own = the facade's structure-specific \c standard keys, written FIRST.  An off gate is
//! spelled "off" (was the 1e30 sentinel).  Inert outside a run, like every cursor write.
void ReportSCFRecipe(const SCFParams& p, const std::string& accelerator, report::json own=report::json::object())
{
    namespace rpt = qchem::report;
    rpt::json standard=std::move(own);
    standard["nMaxIter"]    = long(p.NMaxIter);
    auto gate=[](const Gate& g){return g ? rpt::json(*g) : rpt::json("off");};
    standard["minDrho"]     = gate(p.MinΔρPerElectron);
    standard["minDFD"]      = gate(p.MinΔFD);
    standard["minDE"]       = gate(p.MinΔE);
    standard["minFD"]       = gate(p.MinFD);
    standard["mixer"]       = std::string(p.KerkerG0 > 0.0 ? "Ker" : "Lin") + (p.PulayDepth > 0 ? "+Pul" : "");
    standard["kerkerG0"]    = p.KerkerG0;
    standard["accelerator"] = accelerator;
    standard["smearingkT"]  = p.SmearingkT;

    rpt::json advanced;
    advanced["pulayDepth"]      = p.PulayDepth;
    advanced["pulayStart"]      = p.PulayStart;
    advanced["useMOM"]          = p.UseMOM;
    advanced["momStartIter"]    = p.MOMStartIter;
    advanced["momSmearPenalty"] = p.MOMSmearPenalty;
    advanced["guard"]["holePersistence"] = p.Guard.HolePersistence;
    advanced["guard"]["maxReleases"]     = p.Guard.MaxReleases;
    // the threads this SCF runs on, as the run's scope (ScopedWorkerThreads) or the process set them -- every facade files them here
    advanced["threads"] = {{"openmp", WorkerThreads()}, {"blas", BlasThreads()}};

    rpt::Set("standard", std::move(standard));
    rpt::Set("advanced", std::move(advanced));
}

//! \brief The run's accuracy target eps and WHO CONSUMES it, written ONE LEVEL IN (the caller's `scf` section) and as one console line
//! (V-FACADE-ACCURACY, user ruling 2026-11-03: every facade carries eps, and one whose terms do not read it yet says so -- never
//! silently).  \a consumedBy names the consumers, or says there are none yet.
void ReportAccuracy(const AccuracyTarget& a, const std::string& consumedBy)
{
    std::ostringstream e; e<<a.eps;   // its own stream: the console's formatting state is whatever the caller left
    std::cout<<"[accuracy] eps = "<<e.str()<<" Ha, consumed by: "<<consumedBy<<std::endl;
    report::Set("accuracy", {{"eps", a.eps}, {"consumedBy", consumedBy}});
}

//! What an atom or molecular run says until the eps inventory (doc/Records/EpsilonInventory.md) wires a term of its to the target.
inline constexpr const char* kNoAccuracyConsumerYet = "none yet (the eps inventory)";

} // namespace qchem
