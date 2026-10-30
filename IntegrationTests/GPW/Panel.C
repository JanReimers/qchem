// File: IntegrationTests/GPW/Panel.C  THE CONVERGENCE PANEL's deck runner (doc/SCFRobustnessPlan.md §4, 2026-10-29).
//
// A panel case is a known-hard SCF with an ITERATION BUDGET and a FIXED POINT: a recipe change that helps one case and slows
// another goes red here, before a campaign finds it.  The deck-shaped cases live in decks/panel/ (self-contained: a basis
// file a deck names is vendored beside it and resolved here), and this module runs one in a scratch directory under a
// qchem::IterationBudget (qchem.SCFCap).  It is a helper, not a test: the cases are in their material files.
module;
#include <cstddef>
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

//! The source tree's BasisSetData/ (for a deck that names an in-tree .bsd by file).
std::filesystem::path BasisDataDir() {return std::filesystem::path(BASISSET_DATA_PATH);}

//! Load decks/panel/<name>; a RELATIVE `basis.file` is resolved against decks/panel/ (the vendored basis beside the deck).
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
