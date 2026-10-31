// File: IntegrationTests/GPW/Panel.C  THE CONVERGENCE PANEL's deck runner (doc/SCFRobustnessPlan.md §4, 2026-10-29).
//
// A panel case is a known-hard SCF with an ITERATION BUDGET and a FIXED POINT: a recipe change that helps one case and slows
// another goes red here, before a campaign finds it.  The deck-shaped cases live in decks/panel/ (self-contained: a basis
// file a deck names is vendored beside it and resolved here), and this module runs one in a scratch directory under a
// qchem::IterationBudget (qchem.SCFCap).  It is a helper, not a test: the cases are in their material files.
module;
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>          // getpid: one scratch directory per process
#include <nlohmann/json.hpp>
export module qchem.Tests.Panel;
export import qchem.Deck;     // deck::RunOutcome
export import qchem.SCFCap;   // IterationBudget

export namespace qchem::panel
{

//! decks/panel/ in the source tree (BASISSET_DATA_PATH's sibling).
std::filesystem::path DeckDir()
{
    return std::filesystem::path(BASISSET_DATA_PATH).parent_path()/"decks"/"panel";
}

//! decks/states/ in the source tree: VENDORED converged SCF states (SolidCalculation::SaveState files, a few hundred kB each).
//! An ENERGY anchor tests the Hamiltonian, not the SCF: evaluated at a vendored converged D (an exact resume, <= 3 iterations) it
//! costs one setup instead of a cold convergence (doc/OpenWork.md §3 "Every `_Long` test re-proves its claim at the cheapest
//! resolution").  A state is REFRESHED only by the cold run that owns it, after the anchor is re-judged (doc/Pins.md pin 10).
std::filesystem::path StateDir()
{
    return std::filesystem::path(BASISSET_DATA_PATH).parent_path()/"decks"/"states";
}

//! The vendored state decks/states/<name>.h5.  THROWS if it is absent: a twin without its state is a broken checkout.
std::string VendoredState(const std::string& name)
{
    const std::filesystem::path p=StateDir()/(name+".h5");
    if (!std::filesystem::exists(p)) throw std::runtime_error("panel: no vendored state "+p.string());
    return p.string();
}

//! Where the cold run that OWNS vendored state \a name writes its converged state (the per-process test temp dir), so that
//! refreshing decks/states/<name>.h5 is a copy of a file the run already wrote.  Pair it with PrintRefresh after the run.
std::string RefreshPath(const std::string& name)
{
    return (std::filesystem::temp_directory_path()/(name+"."+std::to_string(::getpid())+".h5")).string();
}

//! Print the one command that refreshes vendored state \a name from the cold run's \a written file.
void PrintRefresh(const std::string& name, const std::string& written)
{
    std::printf("[vendored state] to refresh after re-judging the anchor: cp %s %s\n", written.c_str(),
                (StateDir()/(name+".h5")).string().c_str());
}

//! The source tree's BasisSetData/ (for a deck that names an in-tree .bsd by file).
std::filesystem::path BasisDataDir() {return std::filesystem::path(BASISSET_DATA_PATH);}

//! Load decks/panel/<name>; a RELATIVE `basis.file` is resolved against decks/panel/ (the vendored basis beside the deck), a RELATIVE
//! `state.restartFrom` .h5 against decks/states/ (a vendored saved state).
nlohmann::json LoadDeck(const std::string& name)
{
    const std::filesystem::path p=DeckDir()/name;
    std::ifstream in(p);
    if (!in) throw std::runtime_error("panel: cannot open the deck "+p.string());
    nlohmann::json deck=nlohmann::json::parse(in);
    if (deck.contains("basis") && deck["basis"].contains("file"))
    {
        const std::filesystem::path f=deck["basis"]["file"].get<std::string>();
        if (f.is_relative()) deck["basis"]["file"]=(DeckDir()/f).string();
    }
    if (deck.contains("state") && deck["state"].contains("restartFrom"))
    {   // a relative .h5 names a VENDORED state (decks/states/); a bare stem stays a revision of the run directory
        const std::string r=deck["state"]["restartFrom"].get<std::string>();
        if (r.size()>3 && r.substr(r.size()-3)==".h5" && std::filesystem::path(r).is_relative()) deck["state"]["restartFrom"]=(StateDir()/r).string();
    }
    return deck;
}

//! What a panel deck run returns: the deck's own outcome and the SCF iterations the budget counted (every stage).
struct Outcome
{
    qchem::deck::RunOutcome run;
    std::size_t             iterations=0;
};

//! Run \a deck in a scratch directory (removed afterwards) under an IterationBudget of \a budget named \a what.  Over budget is a
//! failure of the running test (the cap hook); the caller asserts the fixed point.
Outcome Run(const nlohmann::json& deck, std::size_t budget, const char* what)
{
    namespace fs=std::filesystem;
    qchem::deck::RunSpec spec; qchem::deck::FromJson(deck, spec);
    const fs::path dir=fs::temp_directory_path()/("qchem_panel_"+std::to_string(::getpid()));
    fs::remove_all(dir);
    qchem::deck::Provenance pv; pv.codeVersion="panel";
    Outcome o;
    {
        qchem::IterationBudget b(budget, what);
        o.run=qchem::deck::Run(spec, pv, dir);
        o.iterations=b.Used();
    }
    fs::remove_all(dir);
    return o;
}

} // namespace qchem::panel
