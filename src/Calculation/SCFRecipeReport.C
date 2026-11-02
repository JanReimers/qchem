// File: Calculation/SCFRecipeReport.C  The run report's `scf` header -- "what recipe am I running" -- for every facade.
//
// Its own module (V-FACADE-PARITY, 2026-11-02): the molecular and atomic facades each built their own header (two copies of the
// same knobs, the atom's a subset), and the solid facade filed none.  The header is a function of the SCF KNOBS, which are
// structure-neutral; what is structure-specific (the atom's model / spin / electron count) the facade hands in as its own keys.
module;
#include <string>
export module qchem.SCFRecipeReport;
export import qchem.SCFParams;
export import qchem.Reporting;

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

    rpt::Set("standard", std::move(standard));
    rpt::Set("advanced", std::move(advanced));
}

} // namespace qchem
