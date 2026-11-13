// File: Calculation/tests/FacadeMatrix.C  THE FACADE PARITY MATRIX (doc/OOD-SOLID-Cleanup.md V-FACADE-MATRIX, user 2026-11-02:
// "if the OOD/SOLID decoupling was really working well for us there should be almost nothing that works for solids and not for
// molecules (or the other way around)").
//
// ONE parameterised test over {structure kind} x {capability}: every cell exists, so a capability that lands on one facade and not
// the others is a red cell in the sweep, not a GUI report.  A cell that cannot hold YET is an EXPECTED FAILURE naming its tracker
// row (kExpected below) -- and it is asserted to FAIL, so the day the row's work lands the cell goes red and forces the flip (the
// panel's WILL_FAIL rule, CLAUDE.md "Tests", carried to a unit test).  Never delete a cell; flip it.
//
// The fixtures are the cheapest of each kind that exercise a DFT run with every grid: the O atom (Gaussian, LDA), H2O (dzvp, GTH
// LDA: XC + PP quadratures + the auxiliary fit) and Si Γ (SIPP_SR, densityEcut 20) -- seconds each.
#include <gtest/gtest.h>
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <cmath>
#include <exception>
#include <functional>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
import qchem.Calculation;   // Calculation::MakeLevels (the levels cell)
import qchem.StateView;     // the view cell
import qchem.ScalarFunction;
import qchem.MoleculeCalculation;
import qchem.AtomCalculation;
import qchem.SolidCalculation;
import qchem.Deck;
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.StructureData;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData::SIPP_SR
import qchem.SCFParams;
import qchem.SCFIterator;                       // IterationRow
import qchem.SCFAccelerator.Factory;            // SCFAccelerators::Type, SolidAcceleratorOptions
import qchem.Reporting;
import qchem.Types;

using namespace qchem;
using json  = nlohmann::json;   // the deck's
using rjson = report::json;      // the run report's (ordered)

namespace
{
enum class Kind       { Atom, Molecule, Solid };
enum class Capability
{
    PrepareFilesTheSetup,       //!< Prepare files basis.perIrrep + every grids.* entry in its own build run, before any SCF iteration
    ConvergeFilesTheScfHeader,  //!< a converge's run report holds scf.standard (the "what recipe am I running" header)
    AcceleratorKnobsArrive,     //!< every accelerator knob (the Ladder's switchAt included) reaches the accelerator built
    SaveRestartRoundTrips,      //!< a saved state restarts as an EXACT resume (the saved energy, in <= 3 iterations)
    DeckRoundTrips,             //!< a deck naming the structure builds the run, and its json reloads to the same deck
    LevelsAnswerTheState,       //!< the structure-neutral levels (Calculation::MakeLevels): ascending, Σ weight·f = N, at the k they belong to
    ViewMatchesTheLiveState     //!< a saved state, DRAWN (StateView): rho, m, the levels and the orbitals equal the live run's; a moved atom,
                                //!< another basis, another kind and a missing file are REFUSED as Outcomes, never thrown
};

//! THE EXPECTED FAILURES: the row whose work will make the cell hold; nullptr = the cell must hold now.
const char* ExpectedFailure(Kind k, Capability c)
{
    switch (c)
    {
        // A solid's grids section is the raster LADDER's flat bag (grids.densityEcut, grids.ladder, ...), not one entry per grid
        // with its points: the ladder reporting itself as one grid is V-BS-RASTER's.
        case Capability::PrepareFilesTheSetup:      return k==Kind::Solid ? "V-BS-RASTER" : nullptr;
        case Capability::ConvergeFilesTheScfHeader: return nullptr;
        case Capability::AcceleratorKnobsArrive:    return nullptr;
        case Capability::SaveRestartRoundTrips:     return nullptr;   // V-FACADE-STATE landed for the atom and the molecule, 2026-11-09
        case Capability::DeckRoundTrips:            return k==Kind::Solid ? nullptr : "V-FACADE-OPTS";
        case Capability::LevelsAnswerTheState:      return nullptr;
        case Capability::ViewMatchesTheLiveState:   return nullptr;
    }
    return nullptr;
}

//! Every knob, at a NON-default value, so a knob that never arrives cannot pass by matching its default -- read back from the
//! OBJECT THAT USES IT: each accelerator files what it was BUILT with, one level into the context its builder opened (CLAUDE.md
//! "THE RUN REPORT"): the facade opens one `accelerators` row per accelerator it builds, the factory one `rungs` row per Ladder rung.
const rjson kLadder={{"type","Ladder"}, {"eThresh",1e-7}, {"stall",4}, {"floor",1e-9}, {"switchAt",1e-5},
                     {"rungs", rjson::array({ {{"type","DIIS"}, {"nProj",6}, {"engageBelowFD",0.05}, {"svTol",1e-8}},
                                              {{"type","GDM"},  {"trust",0.2}} })}};

//! Every key of \a want is in \a got with the same value (an object recursively, an array element by element); "" when so.
std::string Holds(const rjson& got, const rjson& want, const std::string& path)
{
    if (want.is_object())
    {
        if (!got.is_object()) return path+" is "+got.dump()+", not an object; ";
        std::string why;
        for (auto it=want.begin(); it!=want.end(); ++it)
            why+= got.contains(it.key()) ? Holds(got[it.key()], it.value(), path+"."+it.key()) : path+"."+it.key()+" is absent; ";
        return why;
    }
    if (want.is_array())
    {
        if (!got.is_array() || got.size()!=want.size()) return path+" is "+got.dump()+", asked "+want.dump()+"; ";
        std::string why;
        for (size_t i=0;i<want.size();i++) why+=Holds(got[i], want[i], path+"["+std::to_string(i)+"]");
        return why;
    }
    return got==want ? "" : path+" is "+got.dump()+", asked "+want.dump()+"; ";
}

std::string TempPath(const std::string& leaf) {return (std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_"+leaf)).string();}

//! The newest FILED run whose key starts with \a prefix (e.g. "<ID> build@"); null when none.
const rjson* FiledRun(const std::string& prefix)
{
    const auto& G=report::GlobalReport();
    const rjson* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    return last;
}
std::string FiledKeys()
{
    std::string s;
    for (auto it=report::GlobalReport().begin(); it!=report::GlobalReport().end(); ++it) s+=" '"+it.key()+"'";
    return s.empty() ? " (none)" : s;
}

//! Each verdict below is "" when the capability HOLDS, else WHY it does not.
std::string SetupFiled(const rjson* build, const std::string& prefix, size_t rows)
{
    if (rows>0) return "Prepare ran "+std::to_string(rows)+" SCF iterations";
    if (!build) return "no filed '"+prefix+"' run; filed runs:"+FiledKeys();
    const rjson& b=*build;
    if (!b.contains("basis") || !b["basis"].contains("perIrrep") || b["basis"]["perIrrep"].empty()) return "no basis.perIrrep in the build run";
    for (const rjson& row : b["basis"]["perIrrep"])
        for (const char* k : {"irrep", "lambdaMin", "lambdaMax", "cond"})
            if (!row.contains(k)) return std::string("a basis.perIrrep row has no ")+k;
    if (!b["basis"].contains("conditioning")) return "no basis.conditioning in the build run";
    for (const char* k : {"lambdaMin", "lambdaMinAt", "cond", "condAt", "floor", "wouldDrop"})
        if (!b["basis"]["conditioning"].contains(k)) return std::string("basis.conditioning has no ")+k;
    if (!b.contains("hamiltonian") || b["hamiltonian"].empty()) return "no hamiltonian in the build run";
    const rjson& ham=b["hamiltonian"];
    if (!ham.contains("xcQuadrature")) return "no hamiltonian.xcQuadrature in the build run (a DFT run)";
    // A member of the Hamiltonian's context is a mesh, filed as built, or a term, filing its fit basis.
    std::set<std::string> terms;
    for (auto it=ham.begin(); it!=ham.end(); ++it)
    {
        if (!it->is_object()) return "hamiltonian."+it.key()+" is neither a mesh nor a term";
        if (!it->contains("fitbasis"))
        {
            if (!it->contains("points") || (*it)["points"].get<long>()<=0) return "hamiltonian."+it.key()+" files no points";
            continue;
        }
        terms.insert(it.key());
        for (const auto& fit : (*it)["fitbasis"])
            if (fit.value("role","").empty() || fit.value("name","").empty() || fit.value("nFunctions",0L)<=0)
                return "a hamiltonian."+it.key()+".fitbasis row is incomplete: "+fit.dump();
    }
    for (const char* term : {"Hartree", "Vxc_DFT"}) if (!terms.count(term)) return std::string("no hamiltonian.")+term+".fitbasis";
    if (b.contains("scf") && b["scf"].contains("standard")) return "the build run holds an SCF header";
    return "";
}
std::string ScfHeaderFiled(const rjson* run, const std::string& prefix)
{
    if (!run) return "no filed '"+prefix+"' run; filed runs:"+FiledKeys();
    if (!run->contains("scf") || !(*run)["scf"].contains("standard") || (*run)["scf"]["standard"].empty()) return "the run report has no scf.standard";
    if (!run->contains("result")) return "the run report has no result section";
    for (const char* k : {"converged", "iterations", "charge", "Eelec", "Etot"})
        if (!(*run)["result"].contains(k)) return std::string("result has no ")+k;
    return "";
}
std::string AcceleratorFiled(const rjson* build, const std::string& prefix)
{
    if (!build) return "no filed '"+prefix+"' run; filed runs:"+FiledKeys();
    if (!build->contains("accelerators") || (*build)["accelerators"].empty()) return "the build run files no accelerators row";
    return Holds((*build)["accelerators"].back(), kLadder, "accelerators[last]");
}
std::string DeckBuilds(const json& payload)
{
    try
    {
        deck::RunSpec spec; deck::FromJson(payload, spec);
        const json once=deck::ToJson(spec);
        deck::RunSpec again; deck::FromJson(once, again);
        if (deck::ToJson(again)!=once) return "the deck's json does not reload to the same deck";
        auto calc=deck::Build(spec);
        if (!calc) return "deck::Build returned nothing";
        if (calc->DidConverge()) return "deck::Build ran the SCF";
        return "";
    }
    catch (const std::exception& e) { return std::string("the deck refuses it: ")+e.what(); }
}

//! One structure kind's facade, driven the same way for every capability.
class FacadeProbe
{
public:
    virtual ~FacadeProbe() = default;
    virtual std::string PrepareFilesTheSetup()      = 0;
    virtual std::string ConvergeFilesTheScfHeader() = 0;
    virtual std::string AcceleratorKnobsArrive()    = 0;
    virtual std::string SaveRestartRoundTrips()     = 0;
    virtual std::string DeckRoundTrips()            = 0;
    virtual std::string LevelsAnswerTheState()      = 0;
    virtual std::string ViewMatchesTheLiveState()   = 0;

    std::string Check(Capability c)
    {
        switch (c)
        {
            case Capability::PrepareFilesTheSetup:      return PrepareFilesTheSetup();
            case Capability::ConvergeFilesTheScfHeader: return ConvergeFilesTheScfHeader();
            case Capability::AcceleratorKnobsArrive:    return AcceleratorKnobsArrive();
            case Capability::SaveRestartRoundTrips:     return SaveRestartRoundTrips();
            case Capability::DeckRoundTrips:            return DeckRoundTrips();
            case Capability::LevelsAnswerTheState:      return LevelsAnswerTheState();
            case Capability::ViewMatchesTheLiveState:   return ViewMatchesTheLiveState();
        }
        return "unknown capability";
    }
};

//! The round trip's claim, for every kind (the solid's): an EXACT resume lands on the saved energy to 1e-8 within 3 iterations.
std::string Resumed(bool converged, double E, double Esaved, size_t iters, size_t itersSaved, const std::optional<RestartSummary>& how)
{
    std::ostringstream why;
    if (!converged) why << "the restart did not converge; ";
    if (!how) why << "the restart does not report what it did (RestartedFrom); ";
    else if (!how->exact) { why << "a WARM start, not an exact resume:"; for (const auto& w : how->warm) why << " [" << w << "]"; why << "; "; }
    if (std::fabs(E-Esaved)>1e-8) why << "resumed at " << std::setprecision(12) << E << ", saved " << Esaved << "; ";
    if (iters>3) why << "the resume took " << iters << " iterations (the reference " << itersSaved << ")";
    return why.str();
}

//! The levels' claim, for every kind: a non-empty list, ascending in e, holding the run's \a N electrons (Σ weight·f, to 1e-8), every
//! level at a k the run has (\a gamma: all of them at Γ), and no level holding negative charge.
std::string LevelsHold(const Calculation& c, double N, bool gamma)
{
    const auto levels=c.MakeLevels();
    std::ostringstream why;
    if (levels.empty()) return "no levels";
    double sum=0.0;
    for (size_t i=0;i<levels.size();i++)
    {
        const auto& L=levels[i];
        if (i>0 && L.e<levels[i-1].e) why << "level " << i << " is below level " << i-1 << "; ";
        if (L.occupation<-1e-12) why << "level " << i << " holds " << L.occupation << " electrons; ";
        if (gamma && (L.k.x!=0.0 || L.k.y!=0.0 || L.k.z!=0.0)) why << "level " << i << " is off Γ; ";
        sum+=L.weight*L.occupation;
    }
    if (std::fabs(sum-N)>1e-8) why << "Σ weight·f = " << std::setprecision(12) << sum << ", the run has " << N << " electrons; ";
    return why.str();
}

//! m(r) of an unpolarized live run, for the comparison (the facades answer a zero field or null).
struct ZeroM : public virtual ScalarFunction<double>
{
    double  operator()(const rvec3_t&) const override {return 0.0;}
    rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0.0, 0.0, 0.0);}
};
//! Points around the origin, off every symmetry element (a field that vanished on a mirror would compare as 0 = 0).
const std::vector<rvec3_t>& ProbePoints()
{
    static const std::vector<rvec3_t> pts={{0.0,0.0,0.0}, {0.31,0.17,-0.12}, {1.1,0.53,0.21}, {-0.74,0.41,0.93}, {0.05,-1.2,0.6}, {2.3,-0.4,1.7}};
    return pts;
}
std::string Agree(const char* what, const ScalarFunction<double>& live, const ScalarFunction<double>& drawn)
{
    std::ostringstream why;
    for (const rvec3_t& r : ProbePoints())
    {
        const double a=live(r), b=drawn(r);
        if (!(std::fabs(a-b)<=1e-12*std::max(1.0, std::fabs(a))))
            why << what << " at (" << r.x << "," << r.y << "," << r.z << "): live " << std::setprecision(16) << a << ", drawn " << b << "; ";
    }
    return why.str();
}
//! THE VIEW'S CLAIM, for every kind: a view of the state \a live saved draws what \a live draws -- rho and m at the probe points to
//! 1e-12 (relative), the SAME level list (e, occupation, spin, k, weight, label: bit-equal, it is the saved data), and each orbital the
//! live run samples (\a liveOrbital(i), null = not sampled live) to 1e-12.  Every other level's orbital answers an Outcome, never a throw.
std::string ViewMatches(const StateView& v, const Calculation& live, const ScalarFunction<double>& rho, const ScalarFunction<double>& m,
                        const std::function<const ScalarFunction<double>*(size_t)>& liveOrbital)
{
    std::ostringstream why;
    why << Agree("rho", rho, v.Density()) << Agree("m", m, v.SpinDensity());
    const auto L=live.MakeLevels();
    const auto& V=v.GetLevels();
    if (L.size()!=V.size()) return why.str()+"the view has "+std::to_string(V.size())+" levels, the live run "+std::to_string(L.size());
    for (size_t i=0;i<L.size();i++)
    {
        const auto& a=L[i]; const auto& b=V[i];
        if (a.e!=b.e || a.occupation!=b.occupation || a.spin!=b.spin || a.weight!=b.weight || a.label!=b.label
            || a.k.x!=b.k.x || a.k.y!=b.k.y || a.k.z!=b.k.z)
            { why << "level " << i << ": live " << a.label << " e=" << a.e << " f=" << a.occupation << ", drawn " << b.label << " e=" << b.e
                  << " f=" << b.occupation << "; "; break; }
        try
        {
            auto o=v.Orbital(i);
            if (const ScalarFunction<double>* phi=liveOrbital(i))
            {
                if (!o) why << "level " << i << "'s orbital was refused: " << o.Error().details << "; ";
                else    why << Agree(("orbital "+std::to_string(i)).c_str(), *phi, **o);
            }
        }
        catch (const std::exception& e) { why << "level " << i << "'s orbital THREW: " << e.what() << "; "; }
    }
    try { (void)v.Orbital(L.size()); }   // out of range: an Outcome
    catch (const std::exception& e) { why << "an out-of-range orbital THREW: " << e.what() << "; "; }
    if (!v.Judgement().exact) { why << "the view of the run's own state is not exact:"; for (const auto& w : v.Judgement().warm) why << " [" << w << "]"; why << "; "; }
    return why.str();
}
using ViewCall = std::function<Outcome<std::unique_ptr<StateView>, RestartRefusal>()>;
//! A view that MUST BE REFUSED with \a why (and, when \a mention is not empty, a reason that says it) -- as an Outcome: a throw fails.
std::string Refused(const std::string& what, const ViewCall& view, RestartRefusal::Why why, const std::string& mention="")
{
    try
    {
        auto r=view();
        if (r) return what+" was drawn, not refused; ";
        if (r.Error().why!=why) return what+" was refused for the wrong reason: "+r.Error().details+"; ";
        if (!mention.empty() && r.Error().details.find(mention)==std::string::npos)
            return what+"'s refusal does not say '"+mention+"': "+r.Error().details+"; ";
        return "";
    }
    catch (const std::exception& e) { return what+" THREW instead of refusing: "+e.what()+"; "; }
}
//! The view, opened the way a GUI opens it: inside the caller's run, under its "view" section (the view reports one level in).
Outcome<std::unique_ptr<StateView>, RestartRefusal> Opened(const std::string& run, const ViewCall& view)
{
    report::Run r(run);
    report::Section s("view");
    return view();
}
//! The view's report, filed under the caller's "view" section.
std::string ViewFiled(const std::string& run)
{
    const rjson* f=FiledRun(run+"@");
    if (!f) return "no filed '"+run+"@' run; filed runs:"+FiledKeys();
    if (!f->contains("view")) return "the view filed no 'view' section";
    for (const char* k : {"state", "converged", "energy", "exact", "warm", "electrons", "levels"})
        if (!(*f)["view"].contains(k)) return std::string("the view section has no ")+k;
    return "";
}

class AtomProbe : public FacadeProbe
{
    static constexpr int Z=8;
    static AtomCalcOptions Options() {return {.type=AtomType::Gaussian, .model=Model::LDA};}
    static std::string ID() {return Atom(Z, 0.0, rvec3_t(0,0,0)).ID();}
public:
    std::string PrepareFilesTheSetup() override
    {
        size_t rows=0;
        { report::Run run(ID()+" build"); AtomCalculation::Prepare(Z, 0, Options(), {}, [&](const SCFIterator::IterationRow&){ rows++; }); }   // the probe is the DRIVER: it opens the runs, the facade writes into them (user ruling D)
        return SetupFiled(FiledRun(ID()+" build@"), ID()+" build@", rows);
    }
    std::string ConvergeFilesTheScfHeader() override
    {
        std::unique_ptr<AtomCalculation> calc;
        { report::Run run(ID()); calc=std::make_unique<AtomCalculation>(Z, 0, Options()); }
        if (!calc->IsConverged()) return "the O atom did not converge";
        return ScfHeaderFiled(FiledRun(ID()+"@"), ID()+"@");
    }
    std::string AcceleratorKnobsArrive() override
    {
        // The atom states every knob -- but only through its stringly json escape hatch, in the json door's own key names.
        AtomCalcOptions o=Options();
        o.accelerator={{"type","Ladder"}, {"NProj",6}, {"EMax",0.05}, {"SVTol",1e-8},
                       {"Trust",0.2}, {"ethresh",1e-7}, {"stall",4}, {"floor",1e-9}, {"switchat",1e-5}};
        { report::Run run(ID()+" build"); AtomCalculation::Prepare(Z, 0, o); }
        return AcceleratorFiled(FiledRun(ID()+" build@"), ID()+" build@");
    }
    std::string SaveRestartRoundTrips() override
    {
        const std::string path=TempPath("facade_matrix_o_atom.h5");
        AtomCalculation calc(Z, 0, Options());
        if (!calc.IsConverged()) return "the O atom reference did not converge";
        calc.SaveState(path);
        auto restarted=AtomCalculation::Restart(path, Z, 0, Options());
        std::filesystem::remove(path);
        if (!restarted) return "the restart was refused: "+restarted.Error().details;
        return Resumed((*restarted)->IsConverged(), (*restarted)->Energy(), calc.Energy(), (*restarted)->IterationCount(),
                       calc.IterationCount(), (*restarted)->RestartedFrom());
    }
    std::string DeckRoundTrips()        override {return DeckBuilds({{"structure","B"}});}   // the atom molecules.json carries
    std::string LevelsAnswerTheState() override
    {
        AtomCalculation calc(Z, 0, Options());
        if (!calc.IsConverged()) return "the O atom did not converge";
        return LevelsHold(calc, Z, true);
    }
    //! The atom samples no orbital by level index (its Orbital(i) counts occupied orbitals), so rho, m and the levels are compared,
    //! and every level's orbital must answer.  Refused: another element's spec, a missing file.
    std::string ViewMatchesTheLiveState() override
    {
        const std::string path=TempPath("facade_matrix_view_o_atom.h5");
        AtomCalculation calc(Z, 0, Options());
        if (!calc.IsConverged()) return "the O atom did not converge";
        calc.SaveState(path);
        auto v=Opened("O view", [&]{ return AtomCalculation::ViewState(path, Z, 0, Options()); });
        std::string why;
        if (!v) why="the atom's own state was refused: "+v.Error().details+"; ";
        else
        {
            why+=ViewMatches(**v, calc, calc.Density(), ZeroM(), [](size_t){ return (const ScalarFunction<double>*)nullptr; });
            for (size_t i=0;i<(*v)->GetLevels().size();i++)
                if (!(*v)->Orbital(i)) why+="level "+std::to_string(i)+"'s orbital was refused; ";
            why+=ViewFiled("O view");
        }
        why+=Refused("N's spec on O's state", [&]{ return AtomCalculation::ViewState(path, 7, 0, Options()); }, RestartRefusal::Why::Mismatch);
        std::filesystem::remove(path);
        why+=Refused("a missing file", [&]{ return AtomCalculation::ViewState(path, Z, 0, Options()); }, RestartRefusal::Why::Unreadable);
        return why;
    }
};

class MoleculeProbe : public FacadeProbe
{
    static std::shared_ptr<const MoleculeSetup> Setup(const AcceleratorOptions& acc={})
    {
        return MoleculeSetup::Resolve(StructureData::GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"}, acc);
    }
    static std::string ID() {return StructureData::GetMolecule("H2O").ID();}
public:
    std::string PrepareFilesTheSetup() override
    {
        size_t rows=0;
        { report::Run run(ID()+" build"); MoleculeCalculation::Prepare(Setup(), [&](const SCFIterator::IterationRow&){ rows++; }); }
        return SetupFiled(FiledRun(ID()+" build@"), ID()+" build@", rows);
    }
    std::string ConvergeFilesTheScfHeader() override
    {
        std::unique_ptr<MoleculeCalculation> calc;
        { report::Run run(ID()); calc=std::make_unique<MoleculeCalculation>(Setup()); }
        if (!calc->IsConverged()) return "H2O did not converge";
        return ScfHeaderFiled(FiledRun(ID()+"@"), ID()+"@");
    }
    std::string AcceleratorKnobsArrive() override
    {
        { report::Run run(ID()+" build"); (void)MoleculeCalculation::Prepare(Setup({.nProj=6, .eMax=0.05, .svTol=1e-8, .trust=0.2, .ladderEThresh=1e-7, .ladderStall=4,
                                              .ladderFloor=1e-9, .ladderSwitchAt=1e-5, .type="Ladder"})); }
        return AcceleratorFiled(FiledRun(ID()+" build@"), ID()+" build@");
    }
    std::string SaveRestartRoundTrips() override
    {
        const std::string path=TempPath("facade_matrix_h2o.h5");
        const auto setup=Setup();
        MoleculeCalculation calc(setup);
        auto saved=calc.Result();
        if (!saved) return "the H2O reference did not converge: "+saved.Error().details;
        calc.SaveState(path);
        auto restarted=MoleculeCalculation::Restart(path, Setup());   // a FRESH setup: its mesh is adopted from the state
        std::filesystem::remove(path);
        if (!restarted) return "the restart was refused: "+restarted.Error().details;
        return Resumed((*restarted)->IsConverged(), (*restarted)->Energy(), saved->Energy(), (*restarted)->IterationCount(),
                       saved->IterationCount(), (*restarted)->RestartedFrom());
    }
    std::string DeckRoundTrips()        override {return DeckBuilds({{"structure","H2O"}});}
    //! Plus the molecule's own field-carrying list: index i there is index i here (the binding samples orbital i by it).
    std::string LevelsAnswerTheState() override
    {
        MoleculeCalculation calc(Setup());
        if (!calc.IsConverged()) return "H2O did not converge";
        std::string why=LevelsHold(calc, 8.0, true);   // GTH: O 6 + 2 H
        const auto neutral=calc.MakeLevels();
        const auto& fields=calc.GetLevels();
        if (fields.size()!=neutral.size()) return why+"GetLevels has "+std::to_string(fields.size())+" levels, MakeLevels "+std::to_string(neutral.size());
        for (size_t i=0;i<fields.size();i++)
            if (fields[i].e!=neutral[i].e || fields[i].label!=neutral[i].label || !fields[i].phi)
                return why+"GetLevels["+std::to_string(i)+"] is not MakeLevels["+std::to_string(i)+"]";
        return why;
    }
    //! H2O: rho, m and EVERY orbital (the live GetLevels' fields).  O2 triplet: the spin channels (m is not zero).  Refused: one H
    //! moved by 0.01 bohr (the refusal names the atom and the remedy), the spherical basis, an atom's state, a missing file.
    std::string ViewMatchesTheLiveState() override
    {
        const std::string path=TempPath("facade_matrix_view_h2o.h5"), path2=TempPath("facade_matrix_view_o2.h5"),
                          atomPath=TempPath("facade_matrix_view_atom.h5");
        std::string why;
        {
            MoleculeCalculation calc(Setup());
            if (!calc.IsConverged()) return "H2O did not converge";
            calc.SaveState(path);
            auto v=Opened("H2O view", [&]{ return MoleculeCalculation::ViewState(path, Setup()); });
            if (!v) why+="H2O's own state was refused: "+v.Error().details+"; ";
            else
            {
                const auto& fields=calc.GetLevels();
                why+=ViewMatches(**v, calc, calc.Density(), calc.GetSpinDensity(),
                                 [&](size_t i){ return i<fields.size() ? fields[i].phi : nullptr; });
                why+=ViewFiled("H2O view");
            }
        }
        {
            const auto o2=MoleculeSetup::Resolve(StructureData::GetMolecule("O2"),
                                                 {.basis="dzvp", .model=Model::LDA, .multiplicity=3, .pp="gth-lda-valence"});
            MoleculeCalculation calc(o2);
            if (!calc.IsConverged()) return why+"O2 (triplet) did not converge";
            calc.SaveState(path2);
            auto v=MoleculeCalculation::ViewState(path2, o2);
            if (!v) why+="O2's own state was refused: "+v.Error().details+"; ";
            else
            {
                const auto& fields=calc.GetLevels();
                why+=ViewMatches(**v, calc, calc.Density(), calc.GetSpinDensity(),
                                 [&](size_t i){ return i<fields.size() ? fields[i].phi : nullptr; });
                if (std::fabs((*v)->SpinDensity()(rvec3_t(0.3,0.2,0.1)))<1e-6) why+="O2's drawn m is ~0 near the bond: the channels did not arrive; ";
            }
        }
        {
            Molecule moved;
            size_t a=0;
            const auto h2o=StructureData::GetMolecule("H2O");   // held: Sites() is a VIEW, it must not outlive its structure
            for (const Site& site : h2o.Sites())
            {
                Site m=site;
                if (a++==1) m.R.x+=0.01;
                moved.Insert(m);
            }
            const auto movedSetup=MoleculeSetup::Resolve(moved, {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"});
            why+=Refused("H2O with H 1 moved 0.01 bohr", [&]{ return MoleculeCalculation::ViewState(path, movedSetup); },
                         RestartRefusal::Why::Mismatch, "atom 1 is 0.01 bohr");
            why+=Refused("H2O with H 1 moved 0.01 bohr", [&]{ return MoleculeCalculation::ViewState(path, movedSetup); },
                         RestartRefusal::Why::Mismatch, "reopen the run's own spec, or continue from this run");
        }
        {
            const auto spherical=MoleculeSetup::Resolve(StructureData::GetMolecule("H2O"),
                                                        {.basis="dzvp", .model=Model::LDA, .angular=Angular::Spherical, .pp="gth-lda-valence"});
            why+=Refused("the spherical basis", [&]{ return MoleculeCalculation::ViewState(path, spherical); }, RestartRefusal::Why::Mismatch);
        }
        {
            AtomCalculation atom(8, 0, {.type=AtomType::Gaussian, .model=Model::LDA});
            atom.SaveState(atomPath);
            why+=Refused("an atom's state", [&]{ return MoleculeCalculation::ViewState(atomPath, Setup()); }, RestartRefusal::Why::Kind);
        }
        for (const auto& p : {path, path2, atomPath}) std::filesystem::remove(p);
        why+=Refused("a missing file", [&]{ return MoleculeCalculation::ViewState(path, Setup()); }, RestartRefusal::Why::Unreadable);
        return why;
    }
};

class SolidProbe : public FacadeProbe
{
    static SCFParams Params()
    {
        SCFParams p;
        p.NMaxIter=80; p.MinΔρPerElectron=1e-7; p.MinΔE=1e-10; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;   // off (a Gate of 1e30 would be a LIVE gate)
        return p;
    }
    static std::shared_ptr<const SolidSetup> Setup(SCFAccelerators::Type type=SCFAccelerators::Type::Auto,
                                                   const SCFAccelerators::SolidAcceleratorOptions& acc={})
    {
        const auto mat=Materials::Get("Si_diamond");
        SolidCalcOptions o; o.Nelec=8; o.species=mat.species; o.densityEcut=20.0; o.label="matrix"; o.accelerator=type;
        return SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)),
                                   GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o, acc);
    }
    static std::string ID(const SolidSetup& s) {return s.cell->ID();}
public:
    std::string PrepareFilesTheSetup() override
    {
        const auto setup=Setup();
        size_t rows=0;
        { report::Run run(ID(*setup)+" build"); SolidCalculation::Prepare(setup, {{Params(), SCFAccelerators::Type::Auto}}, [&](const auto&){ rows++; }); }
        return SetupFiled(FiledRun(ID(*setup)+" build@"), ID(*setup)+" build@", rows);
    }
    std::string ConvergeFilesTheScfHeader() override
    {
        const auto setup=Setup();
        std::unique_ptr<SolidCalculation> calc;
        { report::Run run(ID(*setup)); calc=std::make_unique<SolidCalculation>(setup, Params()); }
        if (!calc->DidConverge()) return "Si Γ did not converge";
        return ScfHeaderFiled(FiledRun(ID(*setup)+"@"), ID(*setup)+"@");
    }
    std::string AcceleratorKnobsArrive() override
    {
        SCFAccelerators::SolidAcceleratorOptions acc;
        acc.diis.Nproj=6; acc.diis.EngageBelowFD=0.05; acc.diis.SVTol=1e-8;
        acc.ladderGdm.Trust=0.2; acc.ladderEThresh=1e-7; acc.ladderStall=4; acc.ladderFloor=1e-9; acc.ladderSwitchAt=1e-5;
        const auto setup=Setup(SCFAccelerators::Type::Ladder, acc);
        { report::Run run(ID(*setup)+" build"); SolidCalculation::Prepare(setup, {{Params(), SCFAccelerators::Type::Ladder}}); }
        return AcceleratorFiled(FiledRun(ID(*setup)+" build@"), ID(*setup)+" build@");
    }
    std::string SaveRestartRoundTrips() override
    {
        const std::string path=TempPath("facade_matrix_si.h5");
        const auto setup=Setup();
        SolidCalculation calc(setup, Params());
        auto saved=calc.Result();
        if (!saved) return "the Si Γ reference did not converge: "+saved.Error().details;
        calc.SaveState(path);
        auto restarted=SolidCalculation::Restart(path, setup, Params());
        std::filesystem::remove(path);
        if (!restarted) return "the restart was refused: "+restarted.Error().details;
        auto r=(*restarted)->Result();
        if (!r) return "the restart did not converge: "+r.Error().details;
        std::ostringstream why;
        if (std::fabs(r->Energy()-saved->Energy())>1e-8) why << "resumed at " << r->Energy() << ", saved " << saved->Energy() << "; ";
        if (r->IterationCount()>3) why << "the resume took " << r->IterationCount() << " iterations (the reference " << saved->IterationCount() << ")";
        return why.str();
    }
    std::string DeckRoundTrips() override
    {
        return DeckBuilds({{"structure","Si_diamond"}, {"basis",{{"data","SIPP_SR"}}}, {"solid",{{"densityEcut",20.0}}}});
    }
    std::string LevelsAnswerTheState() override
    {
        SolidCalculation calc(Setup(), Params());
        if (!calc.DidConverge()) return "Si Γ did not converge";
        return LevelsHold(calc, 8.0, true);   // two Si, 4 valence each
    }
    //! Si Γ: rho, m and the levels (the solid samples no orbital live; each level's orbital must answer an Outcome).  THE POSITIONS
    //! ROUND TRIP (controller condition 2, 2026-11-13): a deck whose structure is INLINE with an off-symmetry fraction runs and saves,
    //! its json is written and re-read, and the state is drawn on the re-read spec -- an exact resume of the positions, bit for bit.
    std::string ViewMatchesTheLiveState() override
    {
        const std::string path=TempPath("facade_matrix_view_si.h5"), deckPath=TempPath("facade_matrix_view_si_deck.h5");
        std::string why;
        {
            SolidCalculation calc(Setup(), Params());
            auto r=calc.Result();
            if (!r) return "Si Γ did not converge: "+r.Error().details;
            calc.SaveState(path);
            auto v=Opened("Si view", [&]{ return SolidCalculation::ViewState(path, Setup()); });
            if (!v) why+="Si's own state was refused: "+v.Error().details+"; ";
            else
            {
                const ScalarFunction<double>* m=r->SpinDensity();
                ZeroM zero;
                why+=ViewMatches(**v, calc, r->Density(), m ? *m : zero, [](size_t){ return (const ScalarFunction<double>*)nullptr; });
                why+=ViewFiled("Si view");
            }
            std::filesystem::remove(path);
        }
        try
        {
            const json inlineSi={{"lattice",{{"type","CubicF"},{"a",10.26}}},
                                 {"atoms",{{{"el","Si"},{"frac",{0.0,0.0,0.0}}}, {{"el","Si"},{"frac",{0.25,0.2503333333333333,1.0/3.0-0.0833}}}}}};
            const json payload={{"structure",inlineSi}, {"basis",{{"data","SIPP_SR"}}}, {"solid",{{"densityEcut",20.0}}},
                                {"state",{{"save",deckPath}}}};
            deck::RunSpec spec; deck::FromJson(payload, spec);
            auto calc=deck::Build(spec);
            auto ok=calc->Converge(deck::ScheduleToRun(spec));
            if (!ok) return why+"the inline Si deck did not converge";
            const std::string text=deck::ToJson(spec).dump(2);   // what Run writes as the revision
            deck::RunSpec again; deck::FromJson(json::parse(text), again);
            auto v=deck::View(again, deckPath);
            if (!v) why+="the re-read deck's state was REFUSED (the positions did not round-trip?): "+v.Error().details+"; ";
            else if ((*v)->Judgement().maxDisplacement!=0.0)
                why+="the re-read deck's positions moved "+std::to_string((*v)->Judgement().maxDisplacement)+" bohr; ";
            std::filesystem::remove(deckPath);
        }
        catch (const std::exception& e) { why+=std::string("the inline deck round trip threw: ")+e.what()+"; "; }
        return why;
    }
};

std::unique_ptr<FacadeProbe> Probe(Kind k)
{
    switch (k)
    {
        case Kind::Atom:     return std::make_unique<AtomProbe>();
        case Kind::Molecule: return std::make_unique<MoleculeProbe>();
        case Kind::Solid:    return std::make_unique<SolidProbe>();
    }
    return nullptr;
}

const char* Name(Kind k)       {return k==Kind::Atom ? "Atom" : k==Kind::Molecule ? "Molecule" : "Solid";}
const char* Name(Capability c)
{
    switch (c)
    {
        case Capability::PrepareFilesTheSetup:      return "PrepareFilesTheSetup";
        case Capability::ConvergeFilesTheScfHeader: return "ConvergeFilesTheScfHeader";
        case Capability::AcceleratorKnobsArrive:    return "AcceleratorKnobsArrive";
        case Capability::SaveRestartRoundTrips:     return "SaveRestartRoundTrips";
        case Capability::DeckRoundTrips:            return "DeckRoundTrips";
        case Capability::LevelsAnswerTheState:      return "LevelsAnswerTheState";
        case Capability::ViewMatchesTheLiveState:   return "ViewMatchesTheLiveState";
    }
    return "?";
}

//! One cell's coordinates.  Printed as its test-name suffix, so ctest's discovery (which names a value-parameterised case by its
//! PRINTED value) lists "Facades/FacadeMatrix.Cell/Solid_SaveRestartRoundTrips", not the bytes of two enums.
struct CellId { Kind kind; Capability cap; };
std::string Name(const CellId& c) {return std::string(Name(c.kind))+"_"+Name(c.cap);}
void PrintTo(const CellId& c, std::ostream* os) {*os << Name(c);}

std::vector<CellId> AllCells()
{
    std::vector<CellId> cells;
    for (Kind k : {Kind::Atom, Kind::Molecule, Kind::Solid})
        for (Capability c : {Capability::PrepareFilesTheSetup, Capability::ConvergeFilesTheScfHeader, Capability::AcceleratorKnobsArrive,
                             Capability::SaveRestartRoundTrips, Capability::DeckRoundTrips, Capability::LevelsAnswerTheState,
                             Capability::ViewMatchesTheLiveState})
            cells.push_back({k,c});
    return cells;
}

class FacadeMatrix : public testing::TestWithParam<CellId> {};

} // namespace

TEST_P(FacadeMatrix, Cell)
{
    const auto [kind, cap]=GetParam();
    const std::string why=Probe(kind)->Check(cap);
    if (const char* row=ExpectedFailure(kind, cap))
    {
        EXPECT_FALSE(why.empty()) << Name(GetParam()) << " now HOLDS: " << row << " made it pass -- flip the cell in ExpectedFailure";
        std::cout << "EXPECTED FAILURE (" << row << "): " << Name(GetParam()) << ": " << why << std::endl;
    }
    else
        EXPECT_TRUE(why.empty()) << Name(GetParam()) << ": " << why;
}

INSTANTIATE_TEST_SUITE_P(Facades, FacadeMatrix, testing::ValuesIn(AllCells()),
                         [](const testing::TestParamInfo<CellId>& i) {return Name(i.param);});
