// q2probe -- the Q2 molecular Hartree-route table's driver (lane S 2026-11-25; design note lane-s-design-ij-delta.md, the Q2 spec).
//
//   q2probe <deck.json>
//
// An INSTRUMENT, not a gate (CLAUDE.md: a ladder a human reads is a CLIapps probe): a MOLECULE deck (rundeck drives a cell only so far, OpenWork
// "Deck: Run / Vet / View for a molecule or an atom") through MoleculeSetup::Admit + MoleculeCalculation::Prepare, converged.  (The fixed-D
// audit against exact J joins with V-HAM-COMPOSE step 2, where the exact-J Coulomb term lands.)  One machine-readable line:
//   SETUP  = Admit + Prepare (basis, mesh, fit, integral tables built before the first iteration)
//   SCF    = Converge (every iteration, the mesh pilot's climb and its re-convergence included), with the iteration count
//   peak RSS of the process (getrusage).
// The per-iteration J time is in the run's own ledger (the console's bucket lines), printed by the run itself.
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <nlohmann/json.hpp>
import qchem.Deck;
import qchem.MoleculeCalculation;
import qchem.MoleculeSetup;
import qchem.StructureData;
import qchem.Reporting;

int main(int argc, char** argv)
{
    using namespace qchem;
    if (argc<2) { std::cerr<<"usage: q2probe <deck.json>\n"; return 2; }
    try
    {
        report::Run run("q2probe");
        std::ifstream in(argv[1]);
        const nlohmann::json j=nlohmann::json::parse(in);
        deck::RunSpec spec; deck::FromJson(j, spec);
        using clock=std::chrono::steady_clock;
        const auto t0=clock::now();
        auto setup=MoleculeSetup::Admit(StructureData::GetMolecule(spec.structure), spec.opts, spec.accel);
        if (!setup) { std::cerr<<"q2probe: refused: "<<setup.Error().policy<<"\n"; return 2; }
        auto built=MoleculeCalculation::Prepare(setup.Value());
        if (!built) { std::cerr<<"q2probe: Prepare refused\n"; return 2; }
        auto calc=built.TakeValue();
        const auto t1=clock::now();
        (void)calc->Converge();
        const auto t2=clock::now();
        const bool ok=bool(calc->Result());
        rusage ru{}; getrusage(RUSAGE_SELF, &ru);
        const double setupS=std::chrono::duration<double>(t1-t0).count(), scfS=std::chrono::duration<double>(t2-t1).count();
        std::cout<<std::setprecision(12)<<"Q2ROW label="<<spec.opts.label<<" converged="<<ok<<" E="<<calc->Energy()
                 <<" iters="<<calc->IterationCount()<<" setup_s="<<std::setprecision(4)<<setupS<<" scf_s="<<scfS
                 <<" mesh="<<calc->SizedSetup()->xcMesh.nRadial<<"/"<<calc->SizedSetup()->xcMesh.angularDegree
                 <<" pilotRounds="<<calc->PilotCost().rounds<<" peakRSS_MB="<<ru.ru_maxrss/1024.0;

        std::cout<<std::setprecision(12)<<" Eee="<<calc->EnergyTerms()["Eee"];
        std::cout<<std::endl;
        return ok ? 0 : 1;
    }
    catch (const std::exception& e) { std::cerr<<"q2probe: "<<e.what()<<"\n"; return 2; }
}
