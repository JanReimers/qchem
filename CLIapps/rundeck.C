// rundeck -- run a periodic SCF from an INPUT DECK (D-ENV step 6; doc/Records/EnvKnobInventory.md §9).
//
//   rundeck <deck.json> [--set path=value]... [--out DIR] [--code-version V]
//
// The deck names a structure (a key of materials.json) and states the run; `--set` is the only ad-hoc override (applied before the deck
// is resolved, so it lands in the record).  The run writes `<DIR>/<structure>.rNNN.json` BEFORE it starts: the resolved deck + header +
// provenance, the complete record.  `rundeck <structure>.r003.json` reproduces that run; add a `--set` for the one intentional change.
// DIR defaults to the current directory; use  --out $(scripts/rundir qchem6 <Material>)  to put it where D-RUNDATA keeps runs.
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
import qchem.Deck;
import qchem.Environment;   // RetiredEnvironmentSet / WarnRetiredEnvironment (D-ENV step 6a)

extern char** environ;

namespace {
//! INTERIM (until D-ENV steps 6b/6c delete the remaining hooks): the environment overrides the library still HONOURS.  They are listed
//! in the record so a run is never silently different from its deck.  Resources and diagnostics never change a number and are not listed;
//! the variables 6a retired are IGNORED, so they go under \c ignoredEnvironment instead.
std::vector<std::string> ActiveEnvironment()
{
    std::vector<std::string> out;
    for (char** e=environ; e && *e; ++e)
    {
        const std::string kv(*e); const std::string name=kv.substr(0,kv.find('='));
        if (name.rfind("GPW_",0)!=0 && name.rfind("QCHEM_",0)!=0) continue;
        bool retired=false; for (const auto& r : qchem::RetiredEnvironmentSet()) retired = retired || r.name==name;
        if (retired) continue;
        if (name=="GPW_OMP_THREADS" || name=="QCHEM_OPENMP_THREADS" || name=="QCHEM_BLAS_THREADS" || name=="QCHEM_DIAGNOSTICS") continue;
        out.push_back(kv);
    }
    return out;
}
}

int main(int argc, char** argv)
{
    namespace fs=std::filesystem;
    try
    {
        std::string deckPath, out=".", version="unversioned";
        std::vector<std::string> sets;
        for (int i=1; i<argc; ++i)
        {
            const std::string a=argv[i];
            auto next=[&]{ if (++i>=argc) throw std::runtime_error(a+" needs a value"); return std::string(argv[i]); };
            if      (a=="--set")          sets.push_back(next());
            else if (a=="--out")          out=next();
            else if (a=="--code-version") version=next();
            else if (deckPath.empty() && a.rfind("--",0)!=0) deckPath=a;
            else throw std::runtime_error("unexpected argument '"+a+"'");
        }
        if (deckPath.empty()) { std::cerr<<"usage: rundeck <deck.json> [--set path=value]... [--out DIR] [--code-version V]\n"; return 2; }

        nlohmann::json payload=qchem::deck::LoadDeck(deckPath, version);
        for (const auto& s : sets) qchem::deck::ApplySet(payload, s);
        qchem::deck::RunSpec spec;
        qchem::deck::FromJson(payload, spec);                 // strict: an unknown key (or a --set typo) throws here

        qchem::deck::Provenance prov;
        prov.inputDeck=deckPath; prov.overrides=sets; prov.codeVersion=version; prov.activeEnvironment=ActiveEnvironment();
        { std::ostringstream cl; for (int i=0;i<argc;++i) cl<<(i?" ":"")<<argv[i]; prov.commandLine=cl.str(); }
        qchem::WarnRetiredEnvironment();
        for (const auto& r : qchem::RetiredEnvironmentSet()) prov.ignoredEnvironment.push_back(r.name+"  (use "+r.deckKey+")");
        for (const auto& e : prov.activeEnvironment)
            std::cerr<<"[deck] NOTE: environment override still honoured (retired in D-ENV step 6b/6c): "<<e<<"\n";

        const qchem::deck::RunOutcome r=qchem::deck::Run(spec, prov, out);
        std::cout<<"[rundeck] "<<r.summary<<"\n[rundeck] record: "<<r.revision.string()<<"\n";
        if (r.energy) std::cout.precision(10), std::cout<<"[rundeck] Etot = "<<*r.energy<<" Ha\n";
        return r.converged ? 0 : 1;
    }
    catch (const std::exception& e) { std::cerr<<"rundeck: "<<e.what()<<"\n"; return 2; }
}
