// File: Calculation/Imp/Calculation.C  What every facade shares.
module;
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <functional>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
module qchem.Calculation;
import qchem.Reporting;
import qchem.Policy;             // policy::Refuse: a schedule stage a kind cannot run as stated
import qchem.Hamiltonian.Factory;   // PPFunctionalChecks, ModelToXC (CheckPPFunctional)

namespace qchem
{

void Calculation::ReportResult(const std::string& label, bool converged, size_t iterations, double charge, const EnergyBreakdown& E,
                               double energyTarget, size_t pilotIterations)
{
    namespace rpt = qchem::report;
    const double dq=E.Diagnostic("HartreeChargeError"), dE=E.Diagnostic("HartreeChargeEnergy");
    const bool gridHartree=dq!=0.0;   // only a grid route writes it, and a mesh never integrates a density exactly
    {
        rpt::Section result("result");
        rpt::Set("converged",  converged);
        rpt::Set("iterations", (long)iterations);
        if (pilotIterations>0) rpt::Set("pilotIterations", (long)pilotIterations);
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
             <<"["<<label<<"] "<<(converged ? "CONVERGED" : "NOT converged")<<" iters="<<iterations;
    if (pilotIterations>0) std::cout<<" (pilot "<<pilotIterations<<")";
    std::cout<<" charge="<<charge
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

void Calculation::Checkpoint(const std::string& path, const std::string& label) const
{
    if (path.empty()) return;
    try { SaveState(path); }
    catch (const std::exception& e)
    {
        std::cout << "["<<label<<"] WARNING: state NOT saved to "<<path<<": "<<e.what()<<std::endl;
    }
}

std::string SetupEdit::Describe() const
{
    std::ostringstream os;
    const char* sep="";
    if (eps)   { os<<sep<<"eps -> "<<*eps; sep=", "; }
    if (spin)  { os<<sep<<"spin -> "<<(*spin==SpinGroup::Polarized ? "Polarized" : "UnPolarized"); sep=", "; }
    if (kmesh) { os<<sep<<"kmesh -> "<<kmesh->x<<"x"<<kmesh->y<<"x"<<kmesh->z; }
    return os.str();
}

std::string Calculation::Tag() const
{
    static const std::string unlabelled=RunOptions{}.label;
    const std::string& l=Options().label;
    return l==unlabelled ? KindTag() : l;
}

void Calculation::CheckPPFunctional(const Structure& run) const
{
    const Setup& S=GetSetup();
    if (!S.pp) return;
    policy::PPFittedFunctional(Hamiltonian::PPFunctionalChecks(*S.pp, run, Hamiltonian::ModelToXC(*S.opts.model, S.opts.xalpha)));
}

Calculation::StateIdentity Calculation::MakeStateIdentity() const
{
    // An atom's multiplicity is 0 always (its setup refuses a stated one: the configuration decides its spin); a solid's opts.Nelec and
    // Setup::Ne are the one resolved count (SolidSetup's decision 2).
    const Setup& S=GetSetup();
    return {.Nelec=S.Ne, .spin=StateSpin(), .multiplicity=S.opts.multiplicity, .functional=StateFunctional(), .kT=StateKT()};
}

void Calculation::RefuseSetupEdit(const SetupEdit& edit, size_t index) const
{
    policy::Refuse("StageSetup", "schedule stage "+std::to_string(index)+": this "+StateKind()+" cannot apply the setup edit "+edit.Describe(),
                   "no kind carries its density onto an edited setup yet (doc/OpenWork.md \"Setup-changing schedule stage\"): run the edited "
                   "setup as its own calculation, from a saved state (PrepareRestart judges what carries)");
}

void Calculation::AdmitStageOnSetup(const SCFStage& stage, const RunOptions& opts, SCFAccelerators::Type resolved, const std::string& kind)
{
    using T=SCFAccelerators::Type;
    const bool sameAccel = stage.accelerator==T::Auto ? opts.accelerator==T::Auto || opts.accelerator==resolved
                                                      : stage.accelerator==opts.accelerator || stage.accelerator==resolved;
    if (!sameAccel)
        policy::Refuse("StageAccelerator", std::string("a ")+kind+" schedule stage on accelerator "+NameOf(stage.accelerator)
                       +" (the setup runs "+NameOf(resolved)+")",
                       "a "+kind+" runs every stage on the accelerator its setup admitted (doc/OpenWork.md \"per-stage accelerator on a "
                       "molecule or an atom\"): state it in the setup, or name the setup's in each stage");
    if (stage.recipe!=SCFRecipe::Auto && stage.recipe!=opts.recipe)
        policy::Refuse("StageAccelerator", std::string("a ")+kind+" schedule stage on another recipe than its setup's",
                       "a "+kind+" runs every stage under the recipe its setup states: state it in the setup");
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

const RunOptions& Calculation::Options() const {return GetSetup().opts;}
SpinGroup Calculation::StateSpin() const {return GetSetup().opts.spin;}

//====================================================================================================
//  THE TEMPLATE METHOD (tCalculation<T>): Converge, a schedule, one stage, the outcome
//====================================================================================================
template <class T> Outcome<typename tCalculation<T>::Converged, SCFFailure> tCalculation<T>::Converge(const SCFEdit& edit)
{
    SCFParams params=Options().scf;
    if (edit) edit(params);
    (void)ConvergeStage(params);
    FinishSchedule(params);   // a single Converge is a one-stage schedule: the kind's pilot runs after it too (the solid's rasters)
    return Result();
}

template <class T> Outcome<typename tCalculation<T>::Converged, SCFFailure> tCalculation<T>::Converge(const std::vector<SCFStage>& schedule, const SCFEdit& edit)
{
    if (schedule.empty())
        throw std::invalid_argument("Calculation::Converge: an EMPTY schedule has no meaning -- pass at least one stage, or call Converge()");
    for (size_t s=0; s<schedule.size(); ++s)
    {
        SCFStage stage=schedule[s];
        if (edit) edit(stage.params);
        try
        {
            if (!stage.edit.Empty()) ApplySetupEdit(stage.edit, s);
            EnterStage(stage, s);
        }
        catch (const policy::Refused& r)   // a deck can reach it: an Outcome, its policy row already filed
        {
            SCFFailure f; f.why=SCFFailure::Why::StageRefused; f.details="schedule stage "+std::to_string(s)+" REFUSED (policy "+r.Name()+"): "+r.what();
            return Outcome<Converged, SCFFailure>::Fail(std::move(f));
        }
        (void)ConvergeStage(stage.params);
        LeaveStage(s, schedule.size(), stage.params);
    }
    FinishSchedule(schedule.back().params);
    return Result();
}

template <class T> Outcome<typename tCalculation<T>::Converged, SCFFailure> tCalculation<T>::Result() const
{
    using O=Outcome<Converged, SCFFailure>;
    if (std::optional<SCFFailure> f=MakeFailure()) return O::Fail(std::move(*f));
    return O::Ok(Converged(this));
}

template <class T> Outcome<typename tCalculation<T>::Converged, SCFFailure> tCalculation<T>::ConvergeStage(const SCFParams& params)
{
    if (RunStage(params)) Checkpoint(Options().saveStateTo, Tag());
    return Result();
}

// AN EXACT RESUME STARTS PULAY AT ONCE (OpenWork §3, 2026-10-29; the field-mixer twin of "a restart enters the Ladder at its GDM
// rung").  Plain-Kerker priming exists because a history built far from the fixed point is unstable; an exact resume starts AT it, where
// priming only hands a soft mode five undamped-ish steps: beta-MnO2 A1AFM r015 walked 0.29 Ha away and took 23 iterations back with
// pulayStart 5, 4 with 0.  The history itself is NOT saved (qchem.SolidState: a WARM restart's old residuals extrapolate onto the old
// fixed point); an exact resume rebuilds it in two steps.  It lived in tSCFIterator::Iterate until 2026-11-19 (a policy inside a library).
template <class T> SCFParams tCalculation<T>::ExactResumeMixing(const SCFParams& params, const SCFParams& resolved)
{
    const bool first=itsExactResumePending;
    itsExactResumePending=false;
    if (!first || resolved.PulayDepth<=0 || resolved.PulayStart<=0) return params;
    // A stated depth keeps its stated start through ResolveMixing; with no depth stated the start is Auto's recipe.
    const bool stated = params.PulayDepth>0;
    const std::string decided = "pulayStart "+std::to_string(resolved.PulayStart)+" -> 0 on the first stage of an exact resume"
                              + (stated ? " (overrides the stated pulayStart "+std::to_string(resolved.PulayStart)+")" : "");
    const std::string why = "an exact resume starts AT the fixed point, and plain-Kerker priming is for a start far from it "
                            "(beta-MnO2 A1AFM r015: 23 iterations with pulayStart 5, 4 with 0)";
    std::cout<<"["<<Tag()<<" policy] ExactResumePulay: "<<decided<<"; why: "<<why<<std::endl;
    policy::Report("ExactResumePulay", decided, policy::Verdict::Choose, why, std::nullopt,
                   report::json{{"pulayStart", resolved.PulayStart}, {"stated", stated}});
    SCFParams r=resolved;
    r.PulayStart=0;
    return r;
}

template class tCalculation<double>;   // the atom and the molecule
template class tCalculation<dcmplx>;   // the solid's Bloch blocks

} // namespace qchem
