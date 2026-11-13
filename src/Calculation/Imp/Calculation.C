// File: Calculation/Imp/Calculation.C  What every facade shares.
module;
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
module qchem.Calculation;
import qchem.Reporting;

namespace qchem
{

void Calculation::ReportResult(const std::string& label, bool converged, size_t iterations, double charge, const EnergyBreakdown& E,
                               double energyTarget)
{
    namespace rpt = qchem::report;
    const double dq=E.Diagnostic("HartreeChargeError"), dE=E.Diagnostic("HartreeChargeEnergy");
    const bool gridHartree=dq!=0.0;   // only a grid route writes it, and a mesh never integrates a density exactly
    {
        rpt::Section result("result");
        rpt::Set("converged",  converged);
        rpt::Set("iterations", (long)iterations);
        rpt::Set("charge",     charge);
        rpt::Set("Eelec",      E.GetElectronicEnergy());
        rpt::Set("Etot",       E.GetTotalEnergy());
        if (gridHartree)
        {
            rpt::Set("hartreeChargeError",  dq);
            rpt::Set("hartreeChargeEnergy", dE);
            if (energyTarget>0.0) rpt::Set("hartreeChargeEnergy/target", std::fabs(dE)/energyTarget);
        }
    }
    // Precision STATED, not inherited from whatever the trace left on cout: doc/Benchmark.md compares codes at the 1e-5 Ha level,
    // and anchordiff reads Etot at 1e-10.
    const std::ios_base::fmtflags flags=std::cout.flags();
    const std::streamsize         prec =std::cout.precision();
    std::cout<<std::fixed<<std::setprecision(11)
             <<"["<<label<<"] "<<(converged ? "CONVERGED" : "NOT converged")<<" iters="<<iterations<<" charge="<<charge
             <<" Eelec="<<E.GetElectronicEnergy()<<" Etot="<<std::setprecision(10)<<E.GetTotalEnergy()<<std::setprecision(11)<<"  (";
    const char* sep="";
    for (const auto& [name,shown] : {std::pair{"Kinetic","Ekin"}, {"Een","Een"}, {"Eee","Eee"}, {"Exc","Exc"}, {"Enn","Enn"}, {"E_alphaZ","E_alphaZ"}})
        if (E.Has(name)) { std::cout<<sep<<shown<<"="<<E[name]; sep=" "; }
    std::cout<<")"<<std::endl;
    if (gridHartree)
    {
        std::cout<<std::scientific<<std::setprecision(3)<<"["<<label<<"] Hartree mesh charge error="<<dq<<" energy="<<dE;
        if (energyTarget>0.0) std::cout<<" target="<<energyTarget<<" ratio="<<std::fabs(dE)/energyTarget;
        std::cout<<std::endl;
    }
    std::cout.flags(flags);
    std::cout.precision(prec);
}

//====================================================================================================
//  THE LEVELS: one walk over the state blocks, for every kind
//====================================================================================================
std::vector<Calculation::LevelOrbital> Calculation::WalkLevels() const
{
    std::vector<LevelOrbital> out;
    for (const StateBlock& b : StateBlocks(true))
    {
        if (!b.orbitals) throw std::logic_error("Calculation::WalkLevels: a state block has no orbitals");
        rvec3_t k(0,0,0);
        if (auto it=b.extras.find("k"); it!=b.extras.end() && it->second.size()==3) k=rvec3_t(it->second[0], it->second[1], it->second[2]);
        for (const auto* o : b.orbitals->Iterate())
            out.push_back({{.e=o->GetEigenEnergy(), .occupation=o->GetOccupation(), .spin=b.irrep.ms, .k=k, .weight=b.weight,
                            .label=o->GetLabel()}, o});
    }
    std::stable_sort(out.begin(), out.end(), [](const LevelOrbital& a, const LevelOrbital& b){ return a.level.e < b.level.e; });
    return out;
}

std::vector<Calculation::Level> Calculation::MakeLevels() const
{
    std::vector<Level> out;
    for (auto& l : WalkLevels()) out.push_back(std::move(l.level));
    return out;
}

//====================================================================================================
//  THE SAVED STATE (V-FACADE-STATE): the shared half of every facade's save and restart
//====================================================================================================
void Calculation::RecordState(StateRecord& r, double kT) const
{
    const StateIdentity id=MakeStateIdentity();
    r.blocks=StateBlocks(r.withOrbitals);
    if (r.withOrbitals) r.summary=MakeStateSummary();
    AddStructure(r.fp, StateStructure());
    AddElectrons(r.fp, id.Nelec, id.spin, id.multiplicity, id.functional);
    AddBlocks(r.fp, r.blocks);
    r.fp.warmNum["scf.kT"]={kT>=0.0 ? kT : id.kT};
}

void Calculation::SaveState(const std::string& path) const
{
    StateRecord r;
    r.withOrbitals=true;
    RecordState(r, -1.0);
    {
        qchem::report::Timed timed("state: save (HDF5)");
        WriteState(path, StateKind(), r.fp, r.summary, r.blocks);
    }
    std::cout << "["<<r.summary.label<<"] state saved to "<<path
              << (r.summary.converged ? "" : "  (NOT converged -- a restart will say so)") << std::endl;
}

Outcome<std::vector<std::string>, RestartRefusal> Calculation::JudgeState(const std::string& path) const
{
    using O=Outcome<std::vector<std::string>, RestartRefusal>;
    auto saved=ReadState(path);
    if (!saved) return O::Fail(saved.Error());
    if (saved->kind!=StateKind())
        return O::Fail({RestartRefusal::Why::Kind, "'"+path+"' holds a "+saved->kind+" state; this is a "+StateKind()+" run"});
    StateRecord now;
    RecordState(now, -1.0);
    return CompareFingerprints(saved->fingerprint, now.fp, StateMigration());
}

Outcome<RestartSummary, RestartRefusal> Calculation::JudgeView(const SavedState& saved) const
{
    using O=Outcome<RestartSummary, RestartRefusal>;
    if (saved.kind!=StateKind())
        return O::Fail({RestartRefusal::Why::Kind, "'"+saved.path+"' holds a "+saved.kind+" state; this is a "+StateKind()+" spec"});
    StateRecord now;
    RecordState(now, -1.0);
    const FingerprintMigration migrate=StateMigration();
    auto warm=CompareFingerprints(saved.fingerprint, now.fp, migrate);
    if (!warm) return O::Fail(warm.Error());
    // MOVED ATOMS: refused, never carried.  Name the atom that moved most, by how much, and what to do.
    const StateFingerprint fp=migrate ? migrate(saved.fingerprint, now.fp) : saved.fingerprint;
    const auto sp=fp.warmNum.find("structure.positions");
    const auto np=std::as_const(now.fp).warmNum.find("structure.positions");
    double worst=0.0;
    size_t atom=0;
    if (sp!=fp.warmNum.end() && np!=now.fp.warmNum.end() && sp->second.size()==np->second.size())
        for (size_t a=0;3*a+2<sp->second.size();a++)
        {
            const double dx=sp->second[3*a]-np->second[3*a], dy=sp->second[3*a+1]-np->second[3*a+1], dz=sp->second[3*a+2]-np->second[3*a+2];
            const double d=std::sqrt(dx*dx+dy*dy+dz*dz);
            if (d>worst) {worst=d; atom=a;}
        }
    if (worst>kViewMoveTol)
    {
        std::ostringstream why;
        why<<std::setprecision(6)<<"atom "<<atom<<" is "<<worst<<" bohr from where the saved state has it (the largest move): the spec "
              "was edited since this state was saved; reopen the run's own spec, or continue from this run";
        return O::Fail({RestartRefusal::Why::Mismatch, why.str()});
    }
    return O::Ok(RestartSummary{.path=saved.path, .exact=warm->empty(), .warm=std::move(*warm), .maxDisplacement=worst,
                                .electrons=Electrons(saved, now.blocks)});
}

Outcome<RestartSummary, RestartRefusal> Calculation::JudgeRestart(SavedState& saved, double kT, const std::string& tag) const
{
    using O=Outcome<RestartSummary, RestartRefusal>;
    auto refuse=[&](RestartRefusal r)
    {
        std::cout << tag << "REFUSED ("<<saved.path<<"): "<<r.details<<std::endl;
        return O::Fail(std::move(r));
    };
    if (saved.kind!=StateKind())
        return refuse({RestartRefusal::Why::Kind, "the file holds a "+saved.kind+" state; this is a "+StateKind()+" run"});
    StateRecord now;
    RecordState(now, kT);
    auto warm=CompareFingerprints(saved.fingerprint, now.fp, StateMigration());
    if (!warm) return refuse(warm.Error());

    // MOVED ATOMS (a geometry step, or a state from a nearby geometry): the saved D is a matrix over the OLD centres' functions,
    // carried onto the new ones as AO coefficients (CarryToMovedAtoms, whose declaration holds the measurement: the N-exact
    // carries all started further from the answer).  Its Tr(D S') is a diagnostic, printed below.
    const FingerprintMigration migrate=StateMigration();
    const double moved=MaxDisplacement(migrate ? migrate(saved.fingerprint, now.fp) : saved.fingerprint, now.fp);
    if (moved>1e-10) CarryToMovedAtoms(saved, now.blocks);
    RestartSummary R{.path=saved.path, .exact=warm->empty(), .warm=std::move(*warm), .maxDisplacement=std::max(0.0, moved),
                     .electrons=Electrons(saved, now.blocks)};
    std::ostringstream os;
    os<<tag<<(R.exact ? "EXACT RESUME" : "WARM START")<<" from "<<saved.path<<" (saved by '"<<saved.summary.label<<"', "
      <<(saved.summary.converged ? "converged" : "NOT converged")<<" after "<<saved.summary.iterations
      <<" iterations, E="<<std::setprecision(10)<<saved.summary.energy<<"; "<<saved.summary.recipe<<")";
    for (const std::string& d : R.warm) os<<"\n"<<tag<<"  differs: "<<d;
    if (moved>1e-10)
        os<<"\n"<<tag<<"  the density is carried onto the moved basis as AO coefficients (the functions follow their atoms)";
    // THE SAVED STATE'S OWN RESIDUAL (gap 4, user 2026-10-29): even an exact resume starts where the state was SAVED, not at the
    // fixed point.  A density-MIXED run (no [F,D] measured) stops with D a mixture, one undamped Fock away from self-consistency, so
    // its first iteration re-opens that residual (amplified on a soft system).
    const StateSummary& ss=saved.summary;
    if (std::isfinite(ss.densityChange))
    {
        os<<"\n"<<tag<<"  saved at Δρ "<<std::setprecision(3)<<ss.densityChange;
        if (std::isfinite(ss.commutator)) os<<", [F,D] "<<ss.commutator;
        else os<<" by a density-MIXED run (no [F,D]): D is a mixture, so the first Fock re-opens a residual of about that size -- "
                 "amplified on a soft system (β-MnO₂ r015: 40x), where even an exact resume iterates";
    }
    else os<<"\n"<<tag<<"  saved residual not recorded (a format-1.0 state, before 2026-10-30, or a run that measured none)";
    os<<"\n"<<tag<<"  the restored density holds "<<std::setprecision(12)<<R.electrons<<" electrons on this basis";
    std::cout<<os.str()<<std::endl;
    return O::Ok(std::move(R));
}

} // namespace qchem
