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
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
import qchem.Calculation;
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
    DeckRoundTrips              //!< a deck naming the structure builds the run, and its json reloads to the same deck
};

//! THE EXPECTED FAILURES: the row whose work will make the cell hold; nullptr = the cell must hold now.
const char* ExpectedFailure(Kind k, Capability c)
{
    switch (c)
    {
        case Capability::PrepareFilesTheSetup:      return k==Kind::Solid ? "V-FACADE-PARITY" : nullptr;
        case Capability::ConvergeFilesTheScfHeader: return k==Kind::Solid ? "V-FACADE-PARITY" : nullptr;
        case Capability::AcceleratorKnobsArrive:    return k==Kind::Solid ? "V-FACADE-PARITY" : nullptr;
        case Capability::SaveRestartRoundTrips:     return k==Kind::Solid ? nullptr : "V-FACADE-STATE";
        case Capability::DeckRoundTrips:            return k==Kind::Solid ? nullptr : "V-FACADE-OPTS";
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
    if (!b.contains("grids") || b["grids"].empty()) return "no grids in the build run";
    if (!b["grids"].contains("xcQuadrature")) return "no grids.xcQuadrature in the build run (a DFT run)";
    for (auto it=b["grids"].begin(); it!=b["grids"].end(); ++it)
        if (!it->is_object() || !it->contains("points") || (*it)["points"].get<long>()<=0) return "grids."+it.key()+" files no points";
    if (b.contains("scf") && b["scf"].contains("standard")) return "the build run holds an SCF header";
    return "";
}
std::string ScfHeaderFiled(const rjson* run, const std::string& prefix)
{
    if (!run) return "no filed '"+prefix+"' run; filed runs:"+FiledKeys();
    if (!run->contains("scf") || !(*run)["scf"].contains("standard") || (*run)["scf"]["standard"].empty()) return "the run report has no scf.standard";
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

    std::string Check(Capability c)
    {
        switch (c)
        {
            case Capability::PrepareFilesTheSetup:      return PrepareFilesTheSetup();
            case Capability::ConvergeFilesTheScfHeader: return ConvergeFilesTheScfHeader();
            case Capability::AcceleratorKnobsArrive:    return AcceleratorKnobsArrive();
            case Capability::SaveRestartRoundTrips:     return SaveRestartRoundTrips();
            case Capability::DeckRoundTrips:            return DeckRoundTrips();
        }
        return "unknown capability";
    }
};

class AtomProbe : public FacadeProbe
{
    static constexpr int Z=8;
    static AtomCalcOptions Options() {return {.type=AtomType::Gaussian, .model=Model::LDA};}
    static std::string ID() {return Atom(Z, 0.0, rvec3_t(0,0,0)).ID();}
public:
    std::string PrepareFilesTheSetup() override
    {
        size_t rows=0;
        auto calc=AtomCalculation::Prepare(Z, 0, Options(), {}, [&](const SCFIterator::IterationRow&){ rows++; });
        return SetupFiled(FiledRun(ID()+" build@"), ID()+" build@", rows);
    }
    std::string ConvergeFilesTheScfHeader() override
    {
        AtomCalculation calc(Z, 0, Options());
        if (!calc.IsConverged()) return "the O atom did not converge";
        return ScfHeaderFiled(FiledRun(ID()+"@"), ID()+"@");
    }
    std::string AcceleratorKnobsArrive() override
    {
        // The atom states every knob -- but only through its stringly json escape hatch, in the json door's own key names.
        AtomCalcOptions o=Options();
        o.accelerator={{"type","Ladder"}, {"NProj",6}, {"EMax",0.05}, {"SVTol",1e-8},
                       {"Trust",0.2}, {"ethresh",1e-7}, {"stall",4}, {"floor",1e-9}, {"switchat",1e-5}};
        auto calc=AtomCalculation::Prepare(Z, 0, o);
        return AcceleratorFiled(FiledRun(ID()+" build@"), ID()+" build@");
    }
    std::string SaveRestartRoundTrips() override {return "AtomCalculation has no SaveState / Restart";}
    std::string DeckRoundTrips()        override {return DeckBuilds({{"structure","B"}});}   // the atom molecules.json carries
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
        auto calc=Calculation::Prepare(Setup(), [&](const SCFIterator::IterationRow&){ rows++; });
        return SetupFiled(FiledRun(ID()+" build@"), ID()+" build@", rows);
    }
    std::string ConvergeFilesTheScfHeader() override
    {
        Calculation calc(Setup());
        if (!calc.IsConverged()) return "H2O did not converge";
        return ScfHeaderFiled(FiledRun(ID()+"@"), ID()+"@");
    }
    std::string AcceleratorKnobsArrive() override
    {
        auto calc=Calculation::Prepare(Setup({.nProj=6, .eMax=0.05, .svTol=1e-8, .trust=0.2, .ladderEThresh=1e-7, .ladderStall=4,
                                              .ladderFloor=1e-9, .ladderSwitchAt=1e-5, .type="Ladder"}));
        return AcceleratorFiled(FiledRun(ID()+" build@"), ID()+" build@");
    }
    std::string SaveRestartRoundTrips() override {return "Calculation has no SaveState / Restart";}
    std::string DeckRoundTrips()        override {return DeckBuilds({{"structure","H2O"}});}
};

class SolidProbe : public FacadeProbe
{
    static SCFParams Params()
    {
        SCFParams p;
        p.NMaxIter=80; p.MinΔρPerElectron=1e-7; p.MinΔE=1e-10; p.MinΔFD=1e30; p.MinVirial=1e30; p.MinFD=1e30;
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
        auto calc=SolidCalculation::Prepare(setup, {{Params(), SCFAccelerators::Type::Auto}}, [&](const auto&){ rows++; });
        return SetupFiled(FiledRun(ID(*setup)+" build@"), ID(*setup)+" build@", rows);
    }
    std::string ConvergeFilesTheScfHeader() override
    {
        const auto setup=Setup();
        SolidCalculation calc(setup, Params());
        if (!calc.DidConverge()) return "Si Γ did not converge";
        return ScfHeaderFiled(FiledRun(ID(*setup)+"@"), ID(*setup)+"@");
    }
    std::string AcceleratorKnobsArrive() override
    {
        SCFAccelerators::SolidAcceleratorOptions acc;
        acc.diis.Nproj=6; acc.diis.EngageBelowFD=0.05; acc.diis.SVTol=1e-8;
        acc.ladderGdm.Trust=0.2; acc.ladderEThresh=1e-7; acc.ladderStall=4; acc.ladderFloor=1e-9; acc.ladderSwitchAt=1e-5;
        const auto setup=Setup(SCFAccelerators::Type::Ladder, acc);
        auto calc=SolidCalculation::Prepare(setup, {{Params(), SCFAccelerators::Type::Ladder}});
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
                             Capability::SaveRestartRoundTrips, Capability::DeckRoundTrips})
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
