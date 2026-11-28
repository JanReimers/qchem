// rundeck -- run a periodic SCF from an INPUT DECK (D-ENV step 6; doc/Records/EnvKnobInventory.md §9).
//
//   rundeck <deck.json> [--set path=value]... [--out DIR] [--code-version V] [--check] [--resolved]
//   --check VETS the deck (resolve + postSCF pre-flight + build the lattice and basis) and exits 0 / 2 -- no SCF, no revision.
//
// The deck names a structure (a key of materials.json) and states the run; `--set` is the only ad-hoc override (applied before the deck
// is resolved, so it lands in the record).  The run writes `<DIR>/<structure>.rNNN.json` BEFORE it starts: the resolved deck + header +
// provenance, the complete record.  `rundeck <structure>.r003.json` reproduces that run; add a `--set` for the one intentional change.
// DIR defaults to the current directory; use  --out $(scripts/rundir qchem6 <Material>)  to put it where D-RUNDATA keeps runs.
//
// RUN CONTROL (OpenWork "Run control", user 2026-11-25): a running deck takes commands at each SCF iteration's end --
//   scripts/runctl <pid|DIR> threads N [blas M] | save-and-stop | checkpoint | stop | status      (scripts/runctl list: every live run)
//   SIGTERM / Ctrl-C / SIGUSR1 = save-and-stop (a SECOND SIGTERM or Ctrl-C kills at once).
// A save with no state.save in the deck goes to <DIR>/states/<deck>.stopped.h5, and the verdict line prints the restart.
// EXIT CODES: 0 converged; 1 not converged / refused; 2 a usage or input error; 75 STOPPED by run control (sysexits EX_TEMPFAIL: resumable --
// not 3, which the landing scripts read as "snaprun REFUSED").
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <unistd.h>   // getpid (the run-control banner)
#include <vector>
#include <nlohmann/json.hpp>
import qchem.Deck;
import qchem.FrontEnd.RunControl;   // FrontEndRunControl: runctl, the signals, the heartbeat
import qchem.Environment;   // RetiredEnvironmentSet / WarnRetiredEnvironment (D-ENV step 6a)
import qchem.Reporting;     // report::Run: the driver opens the run the deck is read in

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
        std::string deckPath, out=".", version="unversioned", controlDir; bool check=false, resolved=false;
        std::vector<std::string> sets;
        for (int i=1; i<argc; ++i)
        {
            const std::string a=argv[i];
            auto next=[&]{ if (++i>=argc) throw std::runtime_error(a+" needs a value"); return std::string(argv[i]); };
            if      (a=="--set")          sets.push_back(next());
            else if (a=="--out")          out=next();
            else if (a=="--code-version") version=next();
            else if (a=="--check")        check=true;
            else if (a=="--resolved")     resolved=true;
            else if (a=="--control-dir")  controlDir=next();
            else if (deckPath.empty() && a.rfind("--",0)!=0) deckPath=a;
            else throw std::runtime_error("unexpected argument '"+a+"'");
        }
        if (deckPath.empty())
        {
            std::cerr<<"usage: rundeck <deck.json> [--set path=value]... [--out DIR] [--code-version V] [--check] [--resolved] [--control-dir D]\n"
                       "  while it runs:  scripts/runctl <pid|DIR> threads N [blas M] | save-and-stop | checkpoint | stop | status\n"
                       "                  SIGTERM / Ctrl-C / SIGUSR1 = save-and-stop (a second SIGTERM or Ctrl-C kills)\n"
                       "  exit: 0 converged, 1 not converged or refused, 2 usage/input error, 75 stopped by run control (resumable)\n";
            return 2;
        }

        nlohmann::json payload=qchem::deck::LoadDeck(deckPath, version);
        for (const auto& s : sets) qchem::deck::ApplySet(payload, s);
        qchem::deck::RunSpec spec;
        {   // the driver opens the run the deck is READ in, so a refusal while reading it files its record (its policy rows) too;
            // deck::Run opens the next one, for the resolve, build and SCF
            qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
            const qchem::report::Run read(fs::path(deckPath).stem().string()+" deck");
            qchem::deck::FromJson(payload, spec);             // strict: an unknown key (or a --set typo) throws here
            qchem::report::ClearConsole();
        }

        if (resolved)   // the deck as it will RUN: read, resolved against the structure files, written back (no basis, no SCF)
        {
            (void)qchem::deck::Resolve(spec);
            std::cout<<qchem::deck::ToJson(spec).dump(2)<<"\n";
            return 0;
        }
        if (check)
        {
            const auto v=qchem::deck::Vet(spec, out);
            std::cout<<"[vet] OK  "<<v.summary<<"\n";
            for (const auto& n : v.notes) std::cout<<"[vet] note: "<<n<<"\n";
            return 0;
        }

        qchem::deck::Provenance prov;
        prov.inputDeck=deckPath; prov.overrides=sets; prov.codeVersion=version; prov.activeEnvironment=ActiveEnvironment();
        { std::ostringstream cl; for (int i=0;i<argc;++i) cl<<(i?" ":"")<<argv[i]; prov.commandLine=cl.str(); }
        qchem::WarnRetiredEnvironment();
        for (const auto& r : qchem::RetiredEnvironmentSet()) prov.ignoredEnvironment.push_back(r.name+"  (use "+r.deckKey+")");
        for (const auto& e : prov.activeEnvironment)
            std::cerr<<"[deck] NOTE: environment override still honoured (retired in D-ENV step 6b/6c): "<<e<<"\n";

        // RUN CONTROL: the control directory is the front end's (HOME read here, never in the library); --control-dir overrides it.
        if (controlDir.empty()) { const char* h=std::getenv("HOME"); controlDir=(fs::path(h ? h : ".")/"Code"/"qchem6-runs"/"control").string(); }
        const std::string deckStem=fs::path(deckPath).stem().string();
        auto control=std::make_shared<qchem::frontend::FrontEndRunControl>(
            "rundeck "+deckPath+" -> "+fs::absolute(out).string(), controlDir, (fs::absolute(out)/"states"/(deckStem+".stopped.h5")).string());
        control->InstallSignalHandlers();
        std::cout<<"[rundeck] run control: pid "<<::getpid()<<", commands via scripts/runctl "<<::getpid()<<" ... ("<<control->ControlFile().string()<<")\n";
        { std::error_code ec; fs::create_directories(fs::absolute(out)/"states", ec); }

        const qchem::deck::RunOutcome r=qchem::deck::Run(spec, prov, out, control);
        std::cout<<"[rundeck] "<<r.summary<<"\n[rundeck] record: "<<r.revision.string()<<"\n";
        if (r.energy) std::cout.precision(10), std::cout<<"[rundeck] Etot = "<<*r.energy<<" Ha\n";
        if (r.stopped)
        {
            std::cout<<"[rundeck] STOPPED by run control (exit 75, resumable).  To resume from the saved state:\n"
                     <<"[rundeck]   rundeck "<<deckPath<<" --out "<<out<<" --set state.restartFrom=<the state path above>\n";
            return 75;
        }
        return r.converged ? 0 : 1;
    }
    catch (const std::exception& e) { std::cerr<<"rundeck: "<<e.what()<<"\n"; return 2; }
}
