// File: Calculation/Imp/Deck.C  See the interface.
module;
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <optional>
#include <vector>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <functional>
#include <initializer_list>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>
module qchem.Deck;

import qchem.BasisSet.Lattice.BasisSet;    // CellImages, GPWTolerances
import qchem.Hamiltonian.Factory;           // VxcFit, HubbardManifold
import qchem.SCFAccelerator.Factory;        // SCFAccelerators::Type
import qchem.SCFRecipe;                     // SCFRecipe, AdmitRecipe, FillRecipeKnobs: `scf.recipe` and its admission at --check
import qchem.SCFParams;                     // ResolveMixing: the mixer's admission, run at resolve so --check refuses it too
import qchem.ChargeDensity.Seed;            // SeedStrategy
import qchem.LASolver;                      // qchem::Ortho
import qchem.LASolver.OverlapConditioning;  // ConditioningFloor (MakeBasisShells' default floor)
import qchem.Types;
import qchem.Materials;
import qchem.SolidCalculation;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.BasisFiles;   // BasisFileOfName / BasisFileChecksum: a basis is named by its file stem
import qchem.Reporting;
import qchem.Pseudopotential.Set;           // PseudopotentialSet -- the run's PP, which the basis file's provenance is checked against
import qchem.Outcome;
import qchem.BasisSet.Lattice.BasisSet;   // VetStageTrim, GPWParams
import qchem.BasisSet.Gaussian.Point.ShellTrim;
import qchem.RunPolicy;
import qchem.Policy;                        // policy::Report: the two-run split files its row
import qchem.PeriodicTable;                 // the material's element symbols, for an inlined basis
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94Import;   // SubsetBasis: the elements a deck inlines
import qchem.XCRhoPolicy;                   // RefuseUnroutedGradient: Vet refuses a GGA as the setup will
import qchem.SetupDecisions;                // RequireVxcAdmitted: Vet refuses a pin-30 cell as the setup will
import qchem.SolidSetup;                    // SolidXCMesh: the cell's XC mesh, defaulted as the setup does
import qchem.Mesh.XCPolicy;                 // ResolveXCMesh (a recipe field under cellKind Auto is refused at resolve)
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.StructureData;   // KindOf
import qchem.MoleculeSetup;        // MoleculeSetup::Admit, MoleculeXCMesh: a molecule deck's setup (P0c)
import qchem.MoleculeCalculation;  // MoleculeCalculation::Prepare
import qchem.AtomCalculation;      // AtomCalculation::Prepare, kAtomMeshDefault: an atom deck's (P0c)

namespace qchem::deck
{
namespace
{

//! Reads ONE JSON object into a struct: each \c Get is optional (a missing key keeps the default), and \c Done REJECTS any key
//! no \c Get asked for -- the typo guard.  \a path names the object in the message (`xc.mesh`).
class Fields
{
public:
    Fields(const json& j, std::string path) : itsJ(j), itsPath(std::move(path))
    {
        if (!j.is_object()) throw std::runtime_error("deck: '"+itsPath+"' must be an object, not "+std::string(j.type_name()));
    }
    template <class T> void Get(const char* k, T& v)
    {
        itsKnown.insert(k);
        if (!itsJ.contains(k)) return;
        try { v=itsJ.at(k).get<T>(); }
        catch (const json::exception& e) { throw std::runtime_error("deck: '"+itsPath+"."+k+"': "+e.what()); }
    }
    //! A convergence gate: a number, or the string "off" (\c std::nullopt).  The retired sentinel -- a number >= 1e30, which every
    //! record wrote before 2026-11-02 -- still reads as off, with a notice, so the campaign decks already on disk keep running.
    void GetGate(const char* k, Gate& v)
    {
        itsKnown.insert(k);
        if (!itsJ.contains(k)) return;
        const json& x=itsJ.at(k);
        if (x.is_string() && x.get<std::string>()=="off") { v=std::nullopt; return; }
        if (!x.is_number())
            throw std::runtime_error("deck: '"+itsPath+"."+k+"' = "+x.dump()+" is not a gate: a threshold (a number) or \"off\"");
        const double t=x.get<double>();
        if (t>=1e30)
        {
            std::cout << "[deck] '" << itsPath << "." << k << "' = " << x.dump() << ": the retired 'off' sentinel, read as \"off\" -- "
                         "spell it \"off\"" << std::endl;
            v=std::nullopt; return;
        }
        v=t;
    }
    template <class E> void GetEnum(const char* k, E& v, std::initializer_list<std::pair<const char*,E>> names)
    {
        itsKnown.insert(k);
        if (!itsJ.contains(k)) return;
        const json& x=itsJ.at(k);
        if (x.is_string())
            for (const auto& n : names) if (x.get<std::string>()==n.first) { v=n.second; return; }
        std::string legal; for (const auto& n : names) legal += std::string(legal.empty()?"":", ")+n.first;
        throw std::runtime_error("deck: '"+itsPath+"."+k+"' = "+x.dump()+" is not one of: "+legal);
    }
    const json& At(const char* k) { itsKnown.insert(k); return itsJ.contains(k) ? itsJ.at(k) : kNull; }
    void Done()
    {
        for (auto it=itsJ.begin(); it!=itsJ.end(); ++it)
            if (!itsKnown.count(it.key()) && it.key().rfind("_",0)!=0)      // a key starting with '_' is a comment (the convention of materials.json)
            {
                std::string legal; for (const auto& k : itsKnown) legal += std::string(legal.empty()?"":", ")+k;
                throw std::runtime_error("deck: unknown key '"+itsPath+"."+it.key()+"' (legal: "+legal+")");
            }
    }
private:
    const json& itsJ; std::string itsPath; std::set<std::string> itsKnown; inline static const json kNull=nullptr;
};

template <class E> std::string NameOf(E v, std::initializer_list<std::pair<const char*,E>> names)
{ for (const auto& n : names) if (n.second==v) return n.first; throw std::runtime_error("deck: enum value has no name"); }

using BasisSet::Gaussian::CellImages;
using ChargeDensity::SeedStrategy;
using AccType = SCFAccelerators::Type;
using qcMesh::UnitCellKind; using qcMesh::RadialKind; using qcMesh::AngularKind;

#define NAMES(T,...) std::initializer_list<std::pair<const char*,T>>{__VA_ARGS__}
const auto kImages = NAMES(CellImages,   {"Periodic",CellImages::Periodic}, {"HomeCellOnly",CellImages::HomeCellOnly});
const auto kOrbitalTables = NAMES(policy::XCOrbitalTables, {"auto",policy::XCOrbitalTables::Auto}, {"dense",policy::XCOrbitalTables::Dense},
                                {"images",policy::XCOrbitalTables::Images});
const auto kAngularForm = NAMES(BasisSet::Gaussian::Angular, {"cartesian",BasisSet::Gaussian::Angular::Cartesian}, {"spherical",BasisSet::Gaussian::Angular::Spherical});
const auto kFamily = NAMES(RunOptions::Solid::OrbitalBasis, {"gaussian",RunOptions::Solid::OrbitalBasis::GPW}, {"planeWave",RunOptions::Solid::OrbitalBasis::PlaneWave});
const auto kCoulFit= NAMES(Hamiltonian::CoulombFit, {"Auto",Hamiltonian::CoulombFit::Auto}, {"Auxiliary",Hamiltonian::CoulombFit::Auxiliary}, {"Delta",Hamiltonian::CoulombFit::Delta}, {"Multipole",Hamiltonian::CoulombFit::Multipole});
const auto kVxc    = NAMES(Hamiltonian::VxcFit, {"Auto",Hamiltonian::VxcFit::Auto}, {"PlaneWave",Hamiltonian::VxcFit::PlaneWave}, {"Delta",Hamiltonian::VxcFit::Delta}, {"Auxiliary",Hamiltonian::VxcFit::Auxiliary});
const auto kRho    = NAMES(Hamiltonian::XCρ, {"exact",Hamiltonian::XCρ::Exact}, {"fit",Hamiltonian::XCρ::Fit});
const auto kAcc    = NAMES(AccType, {"DIIS",AccType::DIIS}, {"GDM",AccType::GDM}, {"Ladder",AccType::Ladder}, {"Null",AccType::Null},
                         {"Auto",AccType::Auto});
const auto kRecipe = NAMES(SCFRecipe, {"Auto",SCFRecipe::Auto}, {"Insulator",SCFRecipe::Insulator},
                         {"MagneticInsulator",SCFRecipe::MagneticInsulator}, {"Metal",SCFRecipe::Metal});
const auto kSignal = NAMES(SCFAccelerators::ScheduleSignal, {"Error",SCFAccelerators::ScheduleSignal::Error}, {"EnergyChange",SCFAccelerators::ScheduleSignal::EnergyChange});
const auto kSeed   = NAMES(SeedStrategy, {"Default",SeedStrategy::Default}, {"CoreGuess",SeedStrategy::CoreGuess}, {"Uniform",SeedStrategy::Uniform},
                                         {"SAD",SeedStrategy::SAD}, {"IonicSAD",SeedStrategy::IonicSAD});
const auto kSpin   = NAMES(SpinGroup, {"UnPolarized",SpinGroup::UnPolarized}, {"Polarized",SpinGroup::Polarized});
const auto kEngine = NAMES(BasisSet::Gaussian::Engine, {"MnD",BasisSet::Gaussian::Engine::MnD}, {"LibCint",BasisSet::Gaussian::Engine::LibCint});
const auto kAtomType = NAMES(AtomType, {"Slater",AtomType::Slater}, {"Gaussian",AtomType::Gaussian}, {"BSpline6",AtomType::BSpline6},
                             {"BSpliner6",AtomType::BSpliner6}, {"Gaussian_RKB",AtomType::Gaussian_RKB}, {"Slater_RKB",AtomType::Slater_RKB});
const auto kPool   = NAMES(BasisSetAccuracy, {"Low",BasisSetAccuracy::Low}, {"Medium",BasisSetAccuracy::Medium}, {"High",BasisSetAccuracy::High});
const auto kXC     = NAMES(Hamiltonian::XC, {"SlaterXalpha",Hamiltonian::XC::SlaterXalpha}, {"DiracVWN",Hamiltonian::XC::DiracVWN}, {"LibXC",Hamiltonian::XC::LibXC},
                            {"LibXCGGA",Hamiltonian::XC::LibXCGGA});
const auto kOrtho  = NAMES(Ortho, {"Cholesky",Cholesky}, {"Eigen",Eigen}, {"SVD",SVD}, {"Auto",Auto}, {"CholeskyPivoted",CholeskyPivoted}, {"Orthonormal",Orthonormal});
const auto kCell   = NAMES(UnitCellKind, {"Uniform",UnitCellKind::Uniform}, {"Becke",UnitCellKind::Becke}, {"Auto",UnitCellKind::Auto});
const auto kResolution = NAMES(Response::HubbardResolution, {"manifold",Response::HubbardResolution::Manifold},
                               {"slot",Response::HubbardResolution::Slot}, {"eigen",Response::HubbardResolution::Eigen});
const auto kLoopEstimator = NAMES(SolidCalculation::HubbardLoop::Estimator, {"ACBN0",SolidCalculation::HubbardLoop::Estimator::ACBN0},
                                  {"linearResponse",SolidCalculation::HubbardLoop::Estimator::LinearResponse});
const auto kProjector = NAMES(Hamiltonian::HubbardProjectorKind, {"siteShells",Hamiltonian::HubbardProjectorKind::SiteShells},
                              {"pseudoAtomic",Hamiltonian::HubbardProjectorKind::PseudoAtomic});
const auto kRadial = NAMES(RadialKind, {"MHL",RadialKind::MHL}, {"Log",RadialKind::Log}, {"Linear",RadialKind::Linear}, {"DE2",RadialKind::DE2});
const auto kAngular= NAMES(AngularKind, {"Lebedev",AngularKind::Lebedev}, {"GaussLegendre",AngularKind::GaussLegendre}, {"none",AngularKind::None});
using Measure = SCFParams::Measure;
using BasisSet::Gaussian::BasisSetData;
const auto kMeasure= NAMES(Measure, {"L1PerElectron",Measure::L1PerElectron}, {"MaxDeltaD",Measure::MaxΔD});
const auto kDegenerateFill = NAMES(SCFParams::DegenerateLevelFill, {"Auto",SCFParams::DegenerateLevelFill::Auto},
    {"OrbitalByOrbital",SCFParams::DegenerateLevelFill::OrbitalByOrbital}, {"Shared",SCFParams::DegenerateLevelFill::Shared});
const auto kMixing = NAMES(SCFParams::Mixing, {"Auto",SCFParams::Mixing::Auto}, {"Linear",SCFParams::Mixing::Linear},
                            {"Kerker",SCFParams::Mixing::Kerker}, {"Pulay",SCFParams::Mixing::Pulay});

//! An INLINE basis as a file the reader can open: written once to <temp>/qchem-inline-bases/<FNV-1a of the text>.bsd (content-
//! addressed, so a second deck with the same text reuses it and the file's checksum IS the text's), via a temporary and a rename.
std::string MaterializedBasis(const std::vector<std::string>& lines)
{
    namespace fs=std::filesystem;
    std::string text;
    for (const std::string& l : lines) text+=l+"\n";
    unsigned long long h=1469598103934665603ULL;
    for (unsigned char c : text) { h^=c; h*=1099511628211ULL; }
    char hex[32]; std::snprintf(hex, sizeof hex, "%016llx", h);
    const fs::path dir=fs::temp_directory_path()/"qchem-inline-bases", path=dir/(std::string(hex)+".bsd");
    if (!fs::exists(path))
    {
        fs::create_directories(dir);
        const fs::path tmp=dir/(std::string(hex)+".tmp"+std::to_string(::getpid()));
        { std::ofstream out(tmp); out<<text; if (!out) throw std::runtime_error("deck: cannot write the inline basis to "+tmp.string()); }
        fs::rename(tmp, path);
    }
    return path.string();
}

//! Whether a stated basis.file is one of the SHIPPED BasisSetData files (a bare name, or a path into that directory): those stay a
//! path; any other file is inlined (user ruling 2026-11-10).
bool ShippedBasisFile(const std::string& file)
{
    if (file.find('/')==std::string::npos) return true;
    std::error_code ec;
    return std::filesystem::equivalent(std::filesystem::path(file).parent_path(),
                                       std::filesystem::path(BasisSet::Gaussian::BasisFile("dzvp.bsd")).parent_path(), ec);
}

//! The deck's basis block as the setup's RECIPE (the trim stated, or vetted once by SolidSetup::Resolve).
GaussianBasisRecipe RecipeOf(const RunSpec::Basis& b)
{
    GaussianBasisRecipe r{.file=b.inlineBsd.empty() ? b.file : MaterializedBasis(b.inlineBsd), .autoTrim=b.autoTrim, .angular=b.angular};
    // A REGISTERED name keeps its typed handle; any other stem in the data directory is read as that file.
    if (!BasisSet::Gaussian::RegisteredBasis(b.data, r.data) && r.file.empty()) r.file=BasisSet::Gaussian::BasisFileOfName(b.data);
    for (const auto& t : b.trim) r.trim.shells.push_back({t.Z, t.l, rvec_t(1,t.alpha)});
    return r;
}

} // anon

//=== GPWTolerances ==============================================================================
json ToJson(const BasisSet::Gaussian::GPWTolerances& t)
{
    return {{"vlocEps",t.vlocEps},{"localPPRelCutoff",t.localPPRelCutoff},{"relFieldSharp",t.relFieldSharp},{"mgridEcuts",t.mgridEcuts},
            {"screenEps",t.screenEps},{"fieldSharp",t.fieldSharp},{"relCutoff",t.relCutoff},{"densityEps",t.densityEps},
            {"pairRoutingEps",t.pairRoutingEps},{"rasterX",t.rasterX},{"vlocRasterX",t.vlocRasterX},{"xcRasterX",t.xcRasterX}};
}
void FromJson(const json& j, BasisSet::Gaussian::GPWTolerances& t)
{
    Fields f(j,"tolerances");
    f.Get("vlocEps",t.vlocEps); f.Get("localPPRelCutoff",t.localPPRelCutoff); f.Get("relFieldSharp",t.relFieldSharp);
    f.Get("mgridEcuts",t.mgridEcuts); f.Get("screenEps",t.screenEps); f.Get("fieldSharp",t.fieldSharp);
    f.Get("relCutoff",t.relCutoff); f.Get("densityEps",t.densityEps); f.Get("pairRoutingEps",t.pairRoutingEps);
    f.Get("rasterX",t.rasterX); f.Get("vlocRasterX",t.vlocRasterX); f.Get("xcRasterX",t.xcRasterX);
    f.Done();
}

//=== SCFParams ==================================================================================
//! Is \a g looser than \a t?  An OFF gate is looser than any threshold.
static bool LooserThan(const Gate& g, double t) {return !g || *g>t;}
static json GateJson(const Gate& g) {return g ? json(*g) : json("off");}   //!< a gate as a deck spells it: a number, or "off"
json ToJson(const SCFParams& p)
{
    return {{"NMaxIter",p.NMaxIter},{"minDeltaRhoPerElectron",GateJson(p.MinΔρPerElectron)},{"deltaRhoMeasure",NameOf(p.Δρmeasure,kMeasure)},{"minDeltaFD",GateJson(p.MinΔFD)},
            {"minDeltaE",GateJson(p.MinΔE)},{"minVirial",GateJson(p.MinVirial)},{"minFD",GateJson(p.MinFD)},{"mixer",NameOf(p.Mixer,kMixing)},{"startingRelaxRo",p.StartingRelaxRo},{"maxRelaxRo",p.MaxRelaxRo},{"minRelaxRo",p.MinRelaxRo},{"mergeTol",p.MergeTol},{"degenerateLevelFill",NameOf(p.DegenerateLevels,kDegenerateFill)},
            {"verbose",p.Verbose},{"xcCuspDeficit",p.XCCuspDeficit},{"kerkerG0",p.KerkerG0},{"useMOM",p.UseMOM},{"pulayDepth",p.PulayDepth},
            {"pulayStart",p.PulayStart},{"momStartIter",p.MOMStartIter},{"smearingkT",p.SmearingkT},{"momSmearPenalty",p.MOMSmearPenalty},
            {"stopOnAccelExhausted",p.StopOnAccelExhausted},
            {"momGuard",{{"holePersistence",p.Guard.HolePersistence},{"maxReleases",p.Guard.MaxReleases}}}};
}
void FromJson(const json& j, SCFParams& p)
{
    Fields f(j,"scf");
    if (j.contains("minDeltaRho"))           // RETIRED 2026-10-22 (D-MINDRHO): never silently re-read with a new meaning
        throw std::runtime_error("deck: scf.minDeltaRho is RETIRED -- the gate is now scf.minDeltaRhoPerElectron, the L1 residual per "
                                 "electron on the run's quadrature (the old number meant something different per mixer)");
    f.Get("NMaxIter",p.NMaxIter); f.GetGate("minDeltaRhoPerElectron",p.MinΔρPerElectron); f.GetEnum("deltaRhoMeasure",p.Δρmeasure,kMeasure);
    f.GetGate("minDeltaFD",p.MinΔFD); f.GetGate("minDeltaE",p.MinΔE); f.GetGate("minVirial",p.MinVirial); f.GetGate("minFD",p.MinFD);
    f.GetEnum("mixer",p.Mixer,kMixing); f.Get("startingRelaxRo",p.StartingRelaxRo); f.Get("maxRelaxRo",p.MaxRelaxRo); f.Get("minRelaxRo",p.MinRelaxRo); f.Get("mergeTol",p.MergeTol); f.GetEnum("degenerateLevelFill",p.DegenerateLevels,kDegenerateFill); f.Get("verbose",p.Verbose);
    f.Get("xcCuspDeficit",p.XCCuspDeficit); f.Get("kerkerG0",p.KerkerG0); f.Get("useMOM",p.UseMOM);
    f.Get("pulayDepth",p.PulayDepth); f.Get("pulayStart",p.PulayStart); f.Get("momStartIter",p.MOMStartIter);
    f.Get("smearingkT",p.SmearingkT); f.Get("momSmearPenalty",p.MOMSmearPenalty); f.Get("stopOnAccelExhausted",p.StopOnAccelExhausted);
    if (j.contains("momGuard"))
    {
        Fields g(f.At("momGuard"),"scf.momGuard");
        g.Get("holePersistence",p.Guard.HolePersistence); g.Get("maxReleases",p.Guard.MaxReleases); g.Done();
    }
    else f.At("momGuard");
    f.At("recipe");   // `scf.recipe`, the stage's material regime: not an SCFParams knob -- the RunSpec reader reads it beside them (RecipeOf)
    f.Done();
}
//! A stage's `scf.recipe` (qchem.SCFRecipe; doc/SCFRobustnessPlan.md R4.1 as ruled for B5): Auto when absent, so every older deck reads as
//! it did (decks stay backward compatible, user convention 2026-11-12).  An unknown name throws, listing the legal ones.
static SCFRecipe RecipeOf(const json& scf, const std::string& path)
{
    SCFRecipe r=SCFRecipe::Auto;
    if (!scf.contains("recipe")) return r;
    const json one{{"recipe",scf.at("recipe")}};
    Fields rf(one,path); rf.GetEnum("recipe",r,kRecipe); rf.Done();
    return r;
}
//! A stage's `scf` block as the record writes it: the knobs and the recipe.
static json StageScfJson(const SCFParams& p, SCFRecipe r) {json j=ToJson(p); j["recipe"]=NameOf(r,kRecipe); return j;}

//=== AcceleratorParams ===================================================================================================
static json DiisJson(const SCFAccelerators::DIISParams& d)
{ return {{"nProj",d.Nproj},{"engageBelowFD",d.EngageBelowFD},{"svTol",d.SVTol},{"svTolRel",d.SVTolRel}}; }
static void DiisFromJson(const json& j, SCFAccelerators::DIISParams& d, const std::string& path)
{
    Fields f(j,path); f.Get("nProj",d.Nproj); f.Get("engageBelowFD",d.EngageBelowFD); f.Get("svTol",d.SVTol); f.Get("svTolRel",d.SVTolRel); f.Done();
}
static json GdmJson(const SCFAccelerators::GDMParams& g)
{
    return {{"engageBelowFD",g.EngageBelowFD},{"trust",g.Trust},{"trustBackoff",g.TrustBackoff},{"trustMin",g.TrustMin},{"pcFloor",g.PCFloor},
            {"aufbauTol",g.AufbauTol},{"aufbauMaxRepicks",g.AufbauMaxRepicks}};
}
static void GdmFromJson(const json& j, SCFAccelerators::GDMParams& g, const std::string& path)
{
    Fields f(j,path); f.Get("engageBelowFD",g.EngageBelowFD); f.Get("trust",g.Trust); f.Get("trustBackoff",g.TrustBackoff); f.Get("trustMin",g.TrustMin);
    f.Get("pcFloor",g.PCFloor); f.Get("aufbauTol",g.AufbauTol); f.Get("aufbauMaxRepicks",g.AufbauMaxRepicks); f.Done();
}
json ToJson(const SCFAccelerators::AcceleratorParams& o)
{
    return {{"diis",DiisJson(o.diis)},{"gdm",GdmJson(o.gdm)},
            {"ladder",{{"gdm",GdmJson(o.ladderGdm)},{"eThresh",o.ladderEThresh},{"stall",o.ladderStall},{"floor",o.ladderFloor},
                       {"switchAt",o.ladderSwitchAt},{"signal",NameOf(o.ladderSignal,kSignal)}}}};
}
void FromJson(const json& j, SCFAccelerators::AcceleratorParams& o)
{
    Fields f(j,"acceleratorOptions");
    if (j.contains("diis")) DiisFromJson(f.At("diis"),o.diis,"acceleratorOptions.diis"); else f.At("diis");
    if (j.contains("gdm"))  GdmFromJson (f.At("gdm"), o.gdm, "acceleratorOptions.gdm");  else f.At("gdm");
    if (j.contains("ladder"))
    {
        const json& l=f.At("ladder");
        Fields lf(l,"acceleratorOptions.ladder");
        if (l.contains("gdm")) GdmFromJson(lf.At("gdm"),o.ladderGdm,"acceleratorOptions.ladder.gdm"); else lf.At("gdm");
        lf.Get("eThresh",o.ladderEThresh); lf.Get("stall",o.ladderStall); lf.Get("floor",o.ladderFloor); lf.Get("switchAt",o.ladderSwitchAt);
        lf.GetEnum("signal",o.ladderSignal,kSignal);
        lf.Done();
    }
    else f.At("ladder");
    f.Done();
}
//! The accelerator block as STATED (the inverse of the reader above: an unset knob is not written, so the kind's default decides on reload).
//! A molecule's or an atom's record writes this; the solid's writes its RESOLVED chain (\c AcceleratorParams).
static json StatedJson(const AcceleratorOptions& o)
{
    auto put=[](json& j, const char* k, const auto& v){ if (v) j[k]=*v; };
    auto gdm=[&](const AcceleratorOptions::GDM& g)
    {
        json j=json::object();
        put(j,"engageBelowFD",g.engageBelowFD); put(j,"trust",g.trust); put(j,"trustBackoff",g.trustBackoff); put(j,"trustMin",g.trustMin);
        put(j,"pcFloor",g.pcFloor); put(j,"aufbauTol",g.aufbauTol); put(j,"aufbauMaxRepicks",g.aufbauMaxRepicks);
        return j;
    };
    json j=json::object(), d=json::object(), l=json::object();
    put(d,"nProj",o.diis.nProj); put(d,"engageBelowFD",o.diis.engageBelowFD); put(d,"svTol",o.diis.svTol); put(d,"svTolRel",o.diis.svTolRel);
    if (!d.empty()) j["diis"]=d;
    if (const json g=gdm(o.gdm); !g.empty()) j["gdm"]=g;
    if (const json g=gdm(o.ladder.gdm); !g.empty()) l["gdm"]=g;
    put(l,"eThresh",o.ladder.eThresh); put(l,"stall",o.ladder.stall); put(l,"floor",o.ladder.floor); put(l,"switchAt",o.ladder.switchAt);
    if (o.ladder.signal) l["signal"]=NameOf(*o.ladder.signal,kSignal);
    if (!l.empty()) j["ladder"]=l;
    return j;
}

// The ACCELERATOR BLOCK (V-FACADE-PARITY P0b-1): a deck states only the knobs it means; the rest are the kind's defaults (a periodic
// run: the GPW production recipe), filled by the setup -- and the RECORD writes the resolved, complete chain (ToJson(RunSpec)).
template <class T> static void GetOpt(Fields& f, const json& j, const char* k, std::optional<T>& v)
{ if (j.contains(k)) { T x{}; f.Get(k,x); v=x; } else f.At(k); }
static void GdmOptFromJson(const json& j, AcceleratorOptions::GDM& g, const std::string& path)
{
    Fields f(j,path);
    GetOpt(f,j,"engageBelowFD",g.engageBelowFD); GetOpt(f,j,"trust",g.trust); GetOpt(f,j,"trustBackoff",g.trustBackoff); GetOpt(f,j,"trustMin",g.trustMin);
    GetOpt(f,j,"pcFloor",g.pcFloor); GetOpt(f,j,"aufbauTol",g.aufbauTol); GetOpt(f,j,"aufbauMaxRepicks",g.aufbauMaxRepicks);
    f.Done();
}
void FromJson(const json& j, AcceleratorOptions& o)
{
    Fields f(j,"acceleratorOptions");
    if (j.contains("diis"))
    {
        const json& d=f.At("diis"); Fields df(d,"acceleratorOptions.diis");
        GetOpt(df,d,"nProj",o.diis.nProj); GetOpt(df,d,"engageBelowFD",o.diis.engageBelowFD); GetOpt(df,d,"svTol",o.diis.svTol); GetOpt(df,d,"svTolRel",o.diis.svTolRel);
        df.Done();
    }
    else f.At("diis");
    if (j.contains("gdm")) GdmOptFromJson(f.At("gdm"),o.gdm,"acceleratorOptions.gdm"); else f.At("gdm");
    if (j.contains("ladder"))
    {
        const json& l=f.At("ladder");
        Fields lf(l,"acceleratorOptions.ladder");
        if (l.contains("gdm")) GdmOptFromJson(lf.At("gdm"),o.ladder.gdm,"acceleratorOptions.ladder.gdm"); else lf.At("gdm");
        GetOpt(lf,l,"eThresh",o.ladder.eThresh); GetOpt(lf,l,"stall",o.ladder.stall); GetOpt(lf,l,"floor",o.ladder.floor); GetOpt(lf,l,"switchAt",o.ladder.switchAt);
        if (l.contains("signal")) { SCFAccelerators::ScheduleSignal sg{}; lf.GetEnum("signal",sg,kSignal); o.ladder.signal=sg; } else lf.At("signal");
        lf.Done();
    }
    else f.At("ladder");
    f.Done();
}

//=== MeshParams =================================================================================
json ToJson(const qcMesh::MeshParams& m)
{
    return {{"radial",NameOf(m.radial,kRadial)},{"nRadial",m.nRadial},{"mhlM",m.mhl_m},{"mhlAlpha",m.mhl_alpha},{"logStart",m.logStart},
            {"logStop",m.logStop},{"de2Alpha",m.de2_alpha},{"de2RMin",m.de2_rMin},{"de2RMax",m.de2_rMax},{"angular",NameOf(m.angular,kAngular)},{"angularDegree",m.angularDegree},{"angRot",m.angRot},
            {"beckeOrder",m.beckeOrder},{"nUniform",m.nUniform},{"eCut",m.eCut},{"relCutoff",m.relCutoff},
            {"cellKind",NameOf(qcMesh::UnitCellKind(m.cellKind),kCell)},{"beckeEps",m.beckeEps}};
}
void FromJson(const json& j, qcMesh::MeshParams& m)
{
    Fields f(j,"xcMesh");
    f.GetEnum("radial",m.radial,kRadial); f.Get("nRadial",m.nRadial); f.Get("mhlM",m.mhl_m); f.Get("mhlAlpha",m.mhl_alpha);
    f.Get("logStart",m.logStart); f.Get("logStop",m.logStop);
    f.Get("de2Alpha",m.de2_alpha); f.Get("de2RMin",m.de2_rMin); f.Get("de2RMax",m.de2_rMax); f.GetEnum("angular",m.angular,kAngular);
    f.Get("angularDegree",m.angularDegree); f.Get("angRot",m.angRot); f.Get("beckeOrder",m.beckeOrder); f.Get("nUniform",m.nUniform);
    f.Get("eCut",m.eCut); f.Get("relCutoff",m.relCutoff); { qcMesh::UnitCellKind k=m.cellKind; f.GetEnum("cellKind",k,kCell); m.cellKind=qcMesh::CellKind::RecipeStated(k); } f.Get("beckeEps",m.beckeEps);
    f.Done();
}

//=== RunPolicySpec ==============================================================================
// Only STATED routes are written (an unset optional is "not stated": the umbrella / default decides), so a record reloads to the same
// policy and a later `--set policy.cp2kCompat=true` still means what it says.  The RESOLVED table goes in provenance.policy.
json ToJson(const RunPolicySpec& p)
{
    json j={{"cp2kCompat",p.cp2kCompat}};
    if (p.xcDMMix) j["xcDMMix"]=*p.xcDMMix;
    if (p.xcDMBoost!=1.0) j["xcDMBoost"]=p.xcDMBoost;
    auto put=[&](const char* k, const std::optional<bool>& v){ if (v) j[k]=*v; };
    put("dmLowRank",p.dmLowRank); put("streamFold",p.streamFold); put("mixRhoM",p.mixRhoM); put("xcFromDM",p.xcFromDM);
    put("imposeSymmetry",p.imposeSymmetry); put("beckeXC",p.beckeXC); put("dAwareScreen",p.dAwareScreen); put("hubbardEigen",p.hubbardEigen);
    return j;
}
void FromJson(const json& j, RunPolicySpec& p)
{
    Fields f(j,"policy");
    f.Get("cp2kCompat",p.cp2kCompat);
    if (j.contains("xcDMMix")) { double m=0; f.Get("xcDMMix",m); p.xcDMMix=m; } else f.At("xcDMMix");
    f.Get("xcDMBoost",p.xcDMBoost);
    auto get=[&](const char* k, std::optional<bool>& v){ if (j.contains(k)) { bool b=false; f.Get(k,b); v=b; } else f.At(k); };
    get("dmLowRank",p.dmLowRank); get("streamFold",p.streamFold); get("mixRhoM",p.mixRhoM); get("xcFromDM",p.xcFromDM);
    get("imposeSymmetry",p.imposeSymmetry); get("beckeXC",p.beckeXC); get("dAwareScreen",p.dAwareScreen); get("hubbardEigen",p.hubbardEigen);
    f.Done();
}

//=== RunOptions ===========================================================================
namespace
{
//! The JSON speaks eV (what the literature quotes); RAM is atomic units.  Writing back must give the SAME double on reload, or a record would not
//! reproduce bit-for-bit: pick the shortest-printing eV value whose conversion lands exactly on \a Ha.
constexpr double kHaToEV = 27.211386245988;
double ToEV(double Ha)
{
    const double x=Ha*kHaToEV;
    double best=x; size_t bestLen=std::string::npos;
    double cand[5]={x, std::nextafter(x,1e300), std::nextafter(x,-1e300), std::nextafter(std::nextafter(x,1e300),1e300), std::nextafter(std::nextafter(x,-1e300),-1e300)};
    for (double c : cand)
        if (c/kHaToEV==Ha) { const size_t len=json(c).dump().size(); if (len<bestLen) { best=c; bestLen=len; } }
    return best;
}
std::vector<double> ToEV(const std::vector<double>& v) { std::vector<double> o; for (double x : v) o.push_back(ToEV(x)); return o; }
std::vector<double> FromEV(const std::vector<double>& v) { std::vector<double> o; for (double x : v) o.push_back(x/kHaToEV); return o; }

json ToJson(const Hamiltonian::HubbardManifold& h)
{
    return {{"site",h.site},{"l",h.l},{"U_eV",ToEV(h.U)},{"Uirrep_eV",ToEV(h.Uirrep)},{"alpha_eV",ToEV(h.alpha)},
            {"alphaUp_eV",ToEV(h.alphaUp)},{"alphaDown_eV",ToEV(h.alphaDown)},{"alphaSlot_eV",ToEV(h.alphaSlot)},
            {"orthoAtomic",h.orthoAtomic},{"projector",NameOf(h.projector,kProjector)}};
    // (atomicFunction is not written: the facade derives it from the pseudo-atom, deterministically.)
}
void FromJson(const json& j, Hamiltonian::HubbardManifold& h, const std::string& path)
{
    Fields f(j,path);
    f.Get("site",h.site); f.Get("l",h.l);
    double u=h.U*kHaToEV; f.Get("U_eV",u); h.U=u/kHaToEV;                                   // eV in the file, a.u. in RAM
    std::vector<double> ui=ToEV(h.Uirrep); f.Get("Uirrep_eV",ui); h.Uirrep=FromEV(ui);
    double al=h.alpha*kHaToEV; f.Get("alpha_eV",al); h.alpha=al/kHaToEV;
    double au=h.alphaUp*kHaToEV;   f.Get("alphaUp_eV",au);   h.alphaUp=au/kHaToEV;     // spin-resolved parts, ADDED to alpha_eV
    double ad=h.alphaDown*kHaToEV; f.Get("alphaDown_eV",ad); h.alphaDown=ad/kHaToEV;
    std::vector<double> as=ToEV(h.alphaSlot); f.Get("alphaSlot_eV",as); h.alphaSlot=FromEV(as);   // per U slot, ADDED to alpha_eV
    // RETIRED 2026-10-25 (HubbardProjectorPlan increment 5): the site-shell contraction.  A record written before then
    // carries the inert defaults ("radial": [], "atomicRadial": false), which still load; asking for a contraction does not.
    std::vector<double> radial; bool atomicRadial=false;
    f.Get("radial",radial); f.Get("atomicRadial",atomicRadial);
    if (!radial.empty() || atomicRadial)
        throw std::runtime_error("deck: "+path+": \"radial\" / \"atomicRadial\" are retired (2026-10-25) -- the +U projector is the "
                                 "pseudo-atom orbital by default (\"projector\": \"pseudoAtomic\"); \"siteShells\" names CP2K's LOWDIN manifold");
    f.Get("orthoAtomic",h.orthoAtomic);
    f.GetEnum("projector",h.projector,kProjector);
    f.Done();
}
} // anon

//! \brief THE RUN'S OPTIONS as a schema-2 deck writes them (V-FACADE-PARITY P0, 2026-11-17): the STRUCTURE-NEUTRAL keys at the top
//! (the spin group and multiplicity, the electron count and charge, the pseudopotential set, the model, the seed, the orthogonalisation,
//! the imposed symmetry, the per-site spins, the accelerator, the label, the accuracy target, the XC choices with their mesh), and the
//! periodic-only ones in the \c solid sub-block.  \c ToJson(RunSpec) merges this object into the deck's top level.
json ToJson(const RunOptions& o)
{
    json species=json::array();
    for (const auto& s : o.species) species.push_back({s.first,s.second});
    json hub=json::array();
    for (const auto& h : o.solid.hubbard) hub.push_back(ToJson(h));
    json solid={{"densityEcut",o.solid.densityEcut},{"cutoffFactor",o.solid.cutoffFactor},{"ladderFactor",o.solid.ladderFactor},
            {"images",NameOf(o.solid.images,kImages)},{"kShift",{o.solid.kShift.x,o.solid.kShift.y,o.solid.kShift.z}},
            {"tolerances",ToJson(o.solid.tolerances)},{"policy",ToJson(o.solid.policy)},{"hubbard",hub},
            {"perKFermiLevels",o.solid.perKFermiLevels},{"xcMeshSymmetric",o.solid.xcMeshSymmetric},
            {"xcOrbitalTables",NameOf(o.solid.xcOrbitalTables,kOrbitalTables)},
            {"numBands",o.solid.numBands},
            {"spinsShareFermi",o.solid.spinsShareFermi},{"greyImposition",o.solid.greyImposition},{"momFromSeed",o.solid.momFromSeed},
            {"siteCharges",o.solid.siteCharges}};
    if (o.solid.grid.eta) solid["accuracyCalibration"]={{"vlocLong",o.solid.grid.eta->vlocLong},{"density",o.solid.grid.eta->density}};   // stated only
    if (o.solid.grid.vlocRoutingEta!=BasisSet::Gaussian::GridAccuracy{}.vlocRoutingEta) solid["accuracyCalibration"]["vlocRouting"]=o.solid.grid.vlocRoutingEta;
    if (o.solid.grid.rungGateC!=BasisSet::Gaussian::GridAccuracy{}.rungGateC) solid["accuracyCalibration"]["rungGateC"]=o.solid.grid.rungGateC;
    if (o.solid.grid.relCutoffSafety!=BasisSet::Gaussian::GridAccuracy{}.relCutoffSafety) solid["accuracyCalibration"]["relCutoffSafety"]=o.solid.grid.relCutoffSafety;
    if (o.solid.grid.samplerEta!=BasisSet::Gaussian::GridAccuracy{}.samplerEta) solid["accuracyCalibration"]["samplerEta"]=o.solid.grid.samplerEta;
    if (o.solid.grid.samplerResidualMargin!=BasisSet::Gaussian::GridAccuracy{}.samplerResidualMargin)
        solid["accuracyCalibration"]["samplerResidualMargin"]=o.solid.grid.samplerResidualMargin;
    if (o.solid.firstIterateMargin!=kFirstIterateMargin) solid["firstIterateMargin"]=o.solid.firstIterateMargin;
    if (o.solid.quotientAllowance!=kQuotientAllowance) solid["quotientAllowance"]=o.solid.quotientAllowance;
    if (o.solid.auditFootprint!=kAuditFootprintCalibration)
        solid["auditFootprint"]={{"eta",o.solid.auditFootprint.eta},{"kappa",o.solid.auditFootprint.kappa},{"bandLow",o.solid.auditFootprint.bandLow},
                                 {"bandHigh",o.solid.auditFootprint.bandHigh}};
    if (o.solid.beckeLadder!=qcMesh::BeckeLadderCalibration{})
        solid["beckeLadder"]={{"angularError53",o.solid.beckeLadder.angularError53},{"angularFactor",o.solid.beckeLadder.angularFactor},
                              {"radialError100",o.solid.beckeLadder.radialError100},{"radialFactor",o.solid.beckeLadder.radialFactor},
                              {"radialShare",o.solid.beckeLadder.radialShare},{"angRot",o.solid.beckeLadder.angRot}};
    return {{"Nelec",o.Nelec},{"charge",o.charge},{"multiplicity",o.multiplicity},{"spin",NameOf(o.spin,kSpin)},
            {"pp",o.pp.value_or("gth-lda-valence")},{"species",species},{"model",Hamiltonian::Name(o.model.value_or(Hamiltonian::Model::LDA))},
            {"seed",NameOf(o.seed,kSeed)},{"ortho",NameOf(o.ortho.value_or(qchem::Cholesky),kOrtho)},{"orthoTol",o.orthoTol},
            {"imposeSymmetry",o.imposeSymmetry},{"spins",o.spins},{"accelerator",NameOf(o.accelerator,kAcc)},
            {"label",o.label},{"saveStateTo",o.saveStateTo},
            {"accuracy",o.accuracy.eps},
            {"xc",{{"fit",NameOf(o.vxcFit,kVxc)},{"rho",NameOf(o.xcRho,kRho)},{"mesh",ToJson(o.xcMesh.value_or(qcMesh::XCMeshDefault()))}}},
            {"solid",solid}};
}

static void SolidBlockFromJson(const json& j, RunOptions& o)
{
    Fields f(j,"solid");
    f.Get("densityEcut",o.solid.densityEcut); f.Get("cutoffFactor",o.solid.cutoffFactor); f.Get("ladderFactor",o.solid.ladderFactor);
    if (j.contains("raster"))
        throw std::runtime_error("deck: 'solid.raster' is retired (C2 step 3, doc/GridEpsilonPlan.md §1): the rasters are SIZED from "
                                 "'solid.accuracy' -- drop the key; to give them yourself state \"accuracy\": 0 and tolerances.rasterX "
                                 "(0 = the ball's own raster, the old BallOnly)");
    f.GetEnum("images",o.solid.images,kImages);
    if (j.contains("kShift"))
    {
        const json& k=f.At("kShift");
        if (!k.is_array() || k.size()!=3) throw std::runtime_error("deck: 'solid.kShift' must be [x,y,z], got "+k.dump());
        o.solid.kShift=rvec3_t(k.at(0).get<double>(),k.at(1).get<double>(),k.at(2).get<double>());
    }
    else f.At("kShift");
    if (j.contains("tolerances")) FromJson(f.At("tolerances"),o.solid.tolerances); else f.At("tolerances");
    f.Get("accuracy",o.accuracy.eps);   // the DEPRECATED alias of the top-level `accuracy` (warned and both-refused by ReadRunOptions)
    if (!(o.accuracy.eps>=0.0)) throw std::runtime_error("deck: 'solid.accuracy' (Ha per atom on an energy run, per matrix element on a geometry run) "
                                                         "must be >= 0 (0 = no target: the rasters tolerances.rasterX ... state)");
    if (j.contains("beckeLadder"))   // STATED (expert): where a geometry run's pilot enters the Becke ladder
    {
        Fields c(f.At("beckeLadder"),"solid.beckeLadder");
        qcMesh::BeckeLadderCalibration& b=o.solid.beckeLadder;
        c.Get("angularError53",b.angularError53); c.Get("angularFactor",b.angularFactor); c.Get("radialError100",b.radialError100);
        c.Get("radialFactor",b.radialFactor); c.Get("radialShare",b.radialShare); c.Get("angRot",b.angRot);
        c.Done();
        if (!(b.angularError53>0.0 && b.radialError100>0.0))
            throw std::runtime_error("deck: 'solid.beckeLadder' errors are measured max|dH_xc| values, > 0");
        if (!(b.angularFactor>0.0 && b.angularFactor<1.0 && b.radialFactor>0.0 && b.radialFactor<1.0 && b.radialShare>0.0 && b.radialShare<=1.0))
            throw std::runtime_error("deck: 'solid.beckeLadder' factors are falls per rung, in (0, 1), and radialShare is in (0, 1]");
    }
    f.Get("firstIterateMargin",o.solid.firstIterateMargin);   // STATED (expert): the first-iterate floor's margin, Ha per atom
    if (!(o.solid.firstIterateMargin>=0.0))
        throw std::runtime_error("deck: 'solid.firstIterateMargin' is Ha per atom below the free-atom reference, >= 0");
    f.Get("quotientAllowance",o.solid.quotientAllowance);   // STATED (expert): the quotient range audit's allowance factor
    if (!(o.solid.quotientAllowance>=0.0))
        throw std::runtime_error("deck: 'solid.quotientAllowance' is a relative matrix error, >= 0");
    if (j.contains("auditFootprint"))   // STATED (expert): the audit RSS policy's calibration (absent = its measured record)
    {
        Fields c(f.At("auditFootprint"),"solid.auditFootprint");
        AuditFootprintCalibration& a=o.solid.auditFootprint;
        c.Get("eta",a.eta); c.Get("kappa",a.kappa); c.Get("bandLow",a.bandLow); c.Get("bandHigh",a.bandHigh);
        c.Done();
        if (!(a.eta>0.0 && a.kappa>=0.0 && a.bandLow>0.0 && a.bandLow<=a.bandHigh))
            throw std::runtime_error("deck: 'solid.auditFootprint' needs eta > 0, kappa >= 0 (bytes per raster point) and 0 < bandLow <= bandHigh");
    }
    if (j.contains("accuracyCalibration"))   // STATED (expert); absent = self-calibrated by the run's pilot audit
    {
        // the raster eta is stated by the block (an empty one = rigorous, eta 1); a block holding ONLY the non-raster calibrations
        // (vlocRouting: V_loc long's routing; rungGateC: the ladder's rung gate; relCutoffSafety: the pair rule's stiffness; samplerEta / samplerResidualMargin: the rho-tilde sampler's) leaves the rasters self-calibrated
        const json& cj=j.at("accuracyCalibration");
        bool onlyNonRaster=cj.is_object() && !cj.empty();
        if (cj.is_object()) for (auto it=cj.begin(); it!=cj.end(); ++it)
            if (it.key()!="vlocRouting" && it.key()!="rungGateC" && it.key()!="relCutoffSafety" && it.key()!="samplerEta" && it.key()!="samplerResidualMargin") onlyNonRaster=false;
        const bool etaStated=!onlyNonRaster;
        BasisSet::Gaussian::AccuracyCalibration eta;
        Fields c(f.At("accuracyCalibration"),"solid.accuracyCalibration");
        c.Get("vlocLong",eta.vlocLong); c.Get("density",eta.density);
        c.Get("vlocRouting",o.solid.grid.vlocRoutingEta);
        c.Get("rungGateC",o.solid.grid.rungGateC);
        c.Get("relCutoffSafety",o.solid.grid.relCutoffSafety);
        c.Get("samplerEta",o.solid.grid.samplerEta);
        c.Get("samplerResidualMargin",o.solid.grid.samplerResidualMargin);
        c.Done();
        if (!(o.solid.grid.rungGateC>=0.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.rungGateC' is the C of the ladder's rung gate, >= 0 (0 = no completion rung)");
        if (!(o.solid.grid.relCutoffSafety>0.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.relCutoffSafety' is the stiffness of the pair->level rule, > 0 (default 2)");
        if (!(o.solid.grid.samplerEta>0.0 && o.solid.grid.samplerEta<=1.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.samplerEta' is the rho-tilde NUFFT's delivered/bound, in (0, 1] (1 = rigorous)");
        if (!(o.solid.grid.samplerResidualMargin>0.0 && o.solid.grid.samplerResidualMargin<=1.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.samplerResidualMargin' is the share of the SCF density gate the sampling may use, in (0, 1]");
        if (!(eta.vlocLong>0.0 && eta.vlocLong<=1.0 && eta.density>0.0 && eta.density<=1.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration' factors are eta = delivered/bound, in (0, 1] (1 = rigorous)");
        if (!(o.solid.grid.vlocRoutingEta>0.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.vlocRouting' is V_loc long's routing calibration eta = delivered/bound, > 0 "
                                     "(a stated 'solid.tolerances.vlocEps' wins over it)");
        if (etaStated) o.solid.grid.eta=eta;
    }
    if (j.contains("policy"))     FromJson(f.At("policy"),o.solid.policy);         else f.At("policy");
    if (j.contains("vxcFit"))             // RETIRED 2026-11-10 (V-XC-ONE-TERM step 3): the fit basis is one of the run's XC choices
        throw std::runtime_error("deck: 'solid.vxcFit' is retired (2026-11-10): the Vxc fit basis moved into the run's 'xc' object -- "
                                 "state \"xc\": {\"fit\": \""+(j.at("vxcFit").is_string() ? j.at("vxcFit").get<std::string>() : std::string("Delta"))
                                 +"\"} (beside \"rho\": \"exact\" | \"fit\")");
    if (j.contains("hubbard"))
    {
        o.solid.hubbard.clear();
        size_t i=0;
        for (const auto& h : f.At("hubbard")) { o.solid.hubbard.emplace_back(); FromJson(h,o.solid.hubbard.back(),"solid.hubbard."+std::to_string(i++)); }
    }
    else f.At("hubbard");
    if (j.contains("globalFermi"))           // RETIRED 2026-10-21 (D-FERMI-DEFAULT): never silently re-read with the opposite default
        throw std::runtime_error("deck: solid.globalFermi is RETIRED -- one Fermi level across the k-mesh is now the DEFAULT for a "
                                 "smeared run; delete the key (true), or set solid.perKFermiLevels true for the old per-k default (false)");
    f.Get("perKFermiLevels",o.solid.perKFermiLevels);
    f.Get("xcMeshSymmetric",o.solid.xcMeshSymmetric);   // false = the plain atom-centred mesh: a displacement path needs it (RunOptions)
    f.GetEnum("xcOrbitalTables",o.solid.xcOrbitalTables,kOrbitalTables);   // dense | images | auto: a memory-vs-speed trade-off (RunOptions)
    if (j.contains("forceComplex"))          // RETIRED 2026-11-19 (OpenWork G4, the user: "we can prove rigorously whether it is needed or not")
        throw std::runtime_error("deck: solid.forceComplex is RETIRED -- real vs complex is decided by the setup from k and the planned "
                                 "responses (a postSCF hubbardLinearResponse, linear-response hubbardLoop or hubbardFiniteDifference builds "
                                 "every block complex); delete the key");
    f.Get("numBands",o.solid.numBands);
    f.Get("spinsShareFermi",o.solid.spinsShareFermi); f.Get("greyImposition",o.solid.greyImposition); f.Get("momFromSeed",o.solid.momFromSeed);
    f.Get("siteCharges",o.solid.siteCharges);
    f.Done();
}

//! \brief THE SCHEMA CHECK: a deck is SCHEMA 2 (V-FACADE-PARITY P0, 2026-11-17: the structure-neutral options at the top, only the
//! periodic-only ones in \c solid).  The schema-1 translator that read an old deck's layout is RETIRED (2026-11-18): the GUI's writer
//! writes schema 2 and migrates a saved schema-1 block at load, and the user ruled that keys may move ("we have no customers right now so
//! removing keys is fine as long as we update all our decks accordingly"; a banked \c .rNNN.json record that stops resolving is accepted).
//! So a stated \c "schema": 1, or a key at its schema-1 place (\c solid.multiplicity, ...), is REFUSED by name, naming its schema-2 home
//! and doc/InputDeck.md §3 -- never read with a guessed meaning.  Returns \a j without its \c schema stamp (the readers' field lists do
//! not carry it).
struct MovedKey { const char* old; const char* home; };   // solid.<old> -> <home> ("xc.mesh" = the xc object's mesh)
constexpr MovedKey kSchema1Places[] = {
    {"Nelec","Nelec"}, {"charge","charge"}, {"multiplicity","multiplicity"}, {"pp","pp"}, {"species","species"}, {"model","model"},
    {"seed","seed"}, {"ortho","ortho"}, {"orthoTol","orthoTol"}, {"imposeSymmetry","imposeSymmetry"}, {"siteSpins","spins"},
    {"accelerator","accelerator"}, {"label","label"}, {"saveStateTo","saveStateTo"}, {"xcMesh","xc.mesh"}};

json RequireSchema2(json j)
{
    if (!j.is_object()) return j;
    if (j.contains("schema"))
    {
        if (j.at("schema")==json(1))
            throw std::runtime_error("deck: schema 1 is RETIRED (2026-11-18) -- the structure-neutral options now sit at the top "
                                     "(solid.multiplicity -> multiplicity, solid.siteSpins -> spins, solid.xcMesh -> xc.mesh, ...) and a "
                                     "schema-1 solid.multiplicity >= 1 is stated as \"spin\": \"Polarized\" beside it; see doc/InputDeck.md §3");
        if (j.at("schema")!=json(2))
            throw std::runtime_error("deck: 'schema' is 2 (or absent), got "+j.at("schema").dump()+"; see doc/InputDeck.md §3");
        j.erase("schema");
    }
    if (j.contains("solid") && j.at("solid").is_object())
        for (const MovedKey& m : kSchema1Places)
            if (j.at("solid").contains(m.old))
                throw std::runtime_error(std::string("deck: 'solid.")+m.old+"' is the RETIRED schema-1 place of '"+m.home+"' -- state it as '"
                                         +m.home+"'"+(std::string(m.old)=="multiplicity" ? " (and \"spin\": \"Polarized\" for a polarized run: "
                                                       "schema 1 read multiplicity >= 1 as polarized)" : "")+"; see doc/InputDeck.md §3");
    return j;
}

//! The radial atom's block (\c RunOptions::Atom), every field; \c xc and the \c accelerator escape hatch only when stated.
static json AtomBlockJson(const RunOptions::Atom& a)
{
    json byL=json::array();
    for (const auto& [l,e] : a.exponentsByL) byL.push_back({{"l",l},{"exponents",e}});
    json j={{"type",NameOf(a.type,kAtomType)},{"basisPool",NameOf(a.basisPool,kPool)},{"N",a.N},{"emin",a.emin},{"emax",a.emax},
            {"exponents",a.exponents},{"ltrim",a.ltrim},{"exponentsByL",byL},{"pseudopotential",a.pseudopotential},{"valence",a.valence},
            {"configuration",a.configuration}};
    if (a.xc) j["xc"]={{"kind",NameOf(a.xc->kind,kXC)},{"alpha",a.xc->alpha},{"libxcId",a.xc->libxcId},{"libxcCorrId",a.xc->libxcCorrId}};
    if (!a.accelerator.is_null()) j["accelerator"]=a.accelerator;
    return j;
}
static void AtomBlockFromJson(const json& j, RunOptions::Atom& a)
{
    Fields f(j,"atom");
    f.GetEnum("type",a.type,kAtomType); f.GetEnum("basisPool",a.basisPool,kPool); f.Get("N",a.N); f.Get("emin",a.emin); f.Get("emax",a.emax);
    f.Get("exponents",a.exponents); f.Get("ltrim",a.ltrim); f.Get("pseudopotential",a.pseudopotential); f.Get("valence",a.valence);
    f.Get("configuration",a.configuration);
    if (j.contains("exponentsByL"))
    {
        a.exponentsByL.clear(); size_t i=0;
        for (const auto& e : f.At("exponentsByL"))
        {
            Fields ef(e,"atom.exponentsByL."+std::to_string(i++)); int l=0; std::vector<double> x;
            ef.Get("l",l); ef.Get("exponents",x); ef.Done();
            a.exponentsByL.emplace_back(l,x);
        }
    }
    else f.At("exponentsByL");
    if (j.contains("xc"))
    {
        Fields xf(f.At("xc"),"atom.xc"); Hamiltonian::XCFunctional x;
        xf.GetEnum("kind",x.kind,kXC); xf.Get("alpha",x.alpha); xf.Get("libxcId",x.libxcId); xf.Get("libxcCorrId",x.libxcCorrId); xf.Done();
        a.xc=x;
    }
    else f.At("xc");
    if (j.contains("accelerator")) a.accelerator=f.At("accelerator"); else f.At("accelerator");
    f.Done();
}

//! The run's options from the (translated) deck: the top-level keys through \a f (the caller's typo guard), the \c solid sub-block.
//! The ONE accuracy target and the XC choices are the RUN's (V-FACADE-ACCURACY, V-XC-ONE-TERM): the top-level `accuracy`, and `xc`
//! with its fit basis, the density XC sees and the mesh.  `solid.accuracy` is accuracy's DEPRECATED alias: read, warned; both refused.
void ReadRunOptions(Fields& f, const json& j, RunOptions& o)
{
    f.Get("Nelec",o.Nelec); f.Get("charge",o.charge); f.Get("multiplicity",o.multiplicity); f.GetEnum("spin",o.spin,kSpin); GetOpt(f,j,"pp",o.pp);
    {
        std::string model=o.model ? Hamiltonian::Name(*o.model) : std::string();
        f.Get("model",model);
        if (!model.empty())   // empty: unstated, the kind's default (SolidOptions)
        {
            const auto m=Hamiltonian::ModelFromName(model);
            if (!m) throw std::runtime_error("deck: 'model' \""+model+"\" is not a model name (\"LDA\", \"PADE\", the GGAs \"PBE\", \"PBEsol\", \"BLYP\"; a molecule or an atom also \"HF\", \"Xalpha\", ...)");
            o.model=*m;
        }
    }
    if (j.contains("species"))
    {
        o.species.clear();
        for (const auto& s : f.At("species"))
        {
            if (!s.is_array() || s.size()!=2) throw std::runtime_error("deck: 'species' entries are [\"Element\", valence], got "+s.dump());
            o.species.emplace_back(s.at(0).get<std::string>(), s.at(1).get<int>());
        }
    }
    else f.At("species");
    f.GetEnum("seed",o.seed,kSeed); if (j.contains("ortho")) { qchem::Ortho x{}; f.GetEnum("ortho",x,kOrtho); o.ortho=x; } else f.At("ortho"); f.Get("orthoTol",o.orthoTol);
    f.Get("imposeSymmetry",o.imposeSymmetry); f.Get("spins",o.spins); f.GetEnum("accelerator",o.accelerator,kAcc);
    f.Get("label",o.label); f.Get("saveStateTo",o.saveStateTo);
    f.Get("xalpha",o.xalpha); f.Get("symmetryTol",o.symmetryTol);   // a molecule's (and the atom's xalpha); a solid refuses them (ForeignOption)
    if (j.contains("coulomb"))   // the Hartree route (a molecule's; an atom refuses auxBasis, a solid any of it -- ForeignOption)
    {
        Fields cf(f.At("coulomb"),"coulomb"); cf.GetEnum("fit",o.coulomb.fit,kCoulFit); cf.Get("auxBasis",o.coulomb.auxBasis); cf.Done();
    }
    else f.At("coulomb");
    if (j.contains("atom")) AtomBlockFromJson(f.At("atom"),o.atom); else f.At("atom");
    if (j.contains("solid")) SolidBlockFromJson(f.At("solid"),o); else f.At("solid");
    const bool solidStates = j.contains("solid") && j.at("solid").is_object() && j.at("solid").contains("accuracy");
    if (j.contains("accuracy"))
    {
        if (solidStates)
            throw std::runtime_error("deck: 'accuracy' and 'solid.accuracy' are both stated -- 'solid.accuracy' is the deprecated spelling of the "
                                     "run's ONE accuracy target; keep the top-level 'accuracy' only");
        double eps=o.accuracy.eps;
        f.Get("accuracy",eps);
        if (!(eps>=0.0)) throw std::runtime_error("deck: 'accuracy' (Ha per atom on an energy run, per matrix element on a geometry run) "
                                                  "must be >= 0 (0 = no target: the rasters tolerances.rasterX ... state)");
        o.accuracy.eps=eps;
    }
    else
    {
        f.At("accuracy");
        if (solidStates)
        {
            const std::string w="[deck] WARNING: 'solid.accuracy' is deprecated -- the accuracy target is the run's: state it as the "
                                "top-level 'accuracy' (V-FACADE-ACCURACY)";
            std::fprintf(stderr,"%s\n",w.c_str());
            report::Log(w);
        }
    }
    if (j.contains("xc"))
    {
        Fields xf(f.At("xc"),"xc"); xf.GetEnum("fit",o.vxcFit,kVxc); xf.GetEnum("rho",o.xcRho,kRho);
        if (j.at("xc").contains("mesh")) { qcMesh::MeshParams m=o.xcMesh.value_or(qcMesh::XCMeshDefault()); FromJson(xf.At("mesh"),m); o.xcMesh=m; } else xf.At("mesh");
        xf.Done();
    }
    else f.At("xc");
}

void FromJson(const json& jIn, RunOptions& o)
{
    const json j=RequireSchema2(jIn);
    Fields f(j,"options");
    ReadRunOptions(f,j,o);
    f.Done();
}

//=== RunSpec ====================================================================================
namespace
{
using PostAction=RunSpec::PostAction;
//! The density gate a response REFERENCE converges to (user ruling 2026-10-28; the [F,D] gate replaces it once [F,D] has one home).
constexpr double kResponseReferenceMinΔρ=1e-9;
const char* PostName(PostAction::Kind k)
{
    switch (k)
    {
    case PostAction::Kind::EstimateHubbardU:        return "estimateHubbardU";
    case PostAction::Kind::HubbardLoop:             return "hubbardLoop";
    case PostAction::Kind::IndependentResponse:     return "independentResponse";
    case PostAction::Kind::HubbardLinearResponse:   return "hubbardLinearResponse";
    case PostAction::Kind::HubbardFiniteDifference: return "hubbardFiniteDifference";
    case PostAction::Kind::GridAudit:               return "gridAudit";
    case PostAction::Kind::EggBox:                  return "eggBox";
    case PostAction::Kind::XCMeshAudit:             return "xcMeshAudit";
    case PostAction::Kind::Displace:                return "displace";
    }
    throw std::logic_error("deck: PostAction kind");
}
//! `nq`: an integer n (the cubic n x n x n shorthand, what most decks want) or [n1,n2,n3] like `kmesh` (a supercell gate's 2x1x1).  Written back
//! in the same two forms, so a cubic deck round-trips unchanged.
json QMeshJson(const ivec3_t& q)
{
    if (q.x==q.y && q.y==q.z) return q.x;
    return json::array({q.x,q.y,q.z});
}
void GetQMesh(Fields& f, ivec3_t& nq, const std::string& path)
{
    const json& x=f.At("nq");
    if (x.is_null()) return;
    if (x.is_number_integer()) { const int n=x.get<int>(); nq=ivec3_t(n,n,n); return; }
    if (x.is_array() && x.size()==3 && x[0].is_number_integer() && x[1].is_number_integer() && x[2].is_number_integer())
    { nq=ivec3_t(x[0].get<int>(),x[1].get<int>(),x[2].get<int>()); return; }
    throw std::runtime_error("deck: '"+path+".nq' must be n or [n1,n2,n3], got "+x.dump());
}
json PostJson(const PostAction& a)
{
    json p=json::object();
    switch (a.kind)
    {
    case PostAction::Kind::EstimateHubbardU: break;
    case PostAction::Kind::HubbardLoop:
        p={{"maxOuter",a.maxOuter},{"tolU_eV",a.tolU_eV},{"estimator",NameOf(a.estimator,kLoopEstimator)}};
        if (a.estimator==SolidCalculation::HubbardLoop::Estimator::LinearResponse)
        { p["nq"]=QMeshJson(a.nq); p["perturb"]=a.perturb; p["tol"]=a.tol; p["maxIter"]=a.maxIter; p["restart"]=a.restart; }
        break;
    case PostAction::Kind::IndependentResponse:     p={{"nq",QMeshJson(a.nq)}}; break;
    case PostAction::Kind::HubbardLinearResponse:   p={{"nq",QMeshJson(a.nq)},{"perturb",a.perturb},{"tol",a.tol},{"maxIter",a.maxIter},{"restart",a.restart},
                                                       {"resolve",NameOf(a.resolve,kResolution)},{"qStars",a.qStars},{"littleGroupK",a.littleGroupK}}; break;
    case PostAction::Kind::HubbardFiniteDifference: p={{"perturb",a.perturb},{"alpha_eV",ToEV(a.alpha)}}; break;
    case PostAction::Kind::GridAudit:               p={{"referenceMargin",a.referenceMargin}}; break;
    case PostAction::Kind::EggBox:                  p={{"shift",{a.shift.x,a.shift.y,a.shift.z}}}; break;
    case PostAction::Kind::XCMeshAudit:             p={{"nRadial",a.xcRefNRadial},{"angularDegree",a.xcRefAngularDegree}}; break;
    case PostAction::Kind::Displace:
    {
        json pts=json::array();
        for (const auto& m : a.points) pts.push_back({{"atoms",m.atoms},{"delta",{m.delta.x,m.delta.y,m.delta.z}}});
        p={{"points",pts}}; break;
    }
    }
    return json::object({{PostName(a.kind),p}});
}
PostAction PostFromJson(const json& j, const std::string& path)
{
    if (!j.is_object() || j.size()!=1)
        throw std::runtime_error("deck: '"+path+"' must be an object with exactly ONE key naming the action, e.g. {\"hubbardLoop\":{...}}; got "+j.dump());
    const std::string name=j.begin().key();
    const json& body=j.begin().value();
    PostAction a;
    Fields f(body,path+"."+name);
    if (name=="estimateHubbardU")        a.kind=PostAction::Kind::EstimateHubbardU;
    else if (name=="hubbardLoop")
    {
        a.kind=PostAction::Kind::HubbardLoop; f.Get("maxOuter",a.maxOuter); f.Get("tolU_eV",a.tolU_eV); f.GetEnum("estimator",a.estimator,kLoopEstimator);
        // The response keys only with the estimator that reads them: under ACBN0 they are unknown keys (Done() throws), never silently ignored.
        if (a.estimator==SolidCalculation::HubbardLoop::Estimator::LinearResponse)
        { GetQMesh(f,a.nq,path+"."+name); f.Get("perturb",a.perturb); f.Get("tol",a.tol); f.Get("maxIter",a.maxIter); f.Get("restart",a.restart); }
    }
    else if (name=="independentResponse"){ a.kind=PostAction::Kind::IndependentResponse; GetQMesh(f,a.nq,path+"."+name); }
    else if (name=="hubbardLinearResponse")
    { a.kind=PostAction::Kind::HubbardLinearResponse; GetQMesh(f,a.nq,path+"."+name); f.Get("perturb",a.perturb); f.Get("tol",a.tol); f.Get("maxIter",a.maxIter); f.Get("restart",a.restart);
      f.GetEnum("resolve",a.resolve,kResolution); f.Get("qStars",a.qStars); f.Get("littleGroupK",a.littleGroupK); }
    else if (name=="hubbardFiniteDifference")
    { a.kind=PostAction::Kind::HubbardFiniteDifference; f.Get("perturb",a.perturb); double eV=ToEV(a.alpha); f.Get("alpha_eV",eV); a.alpha=eV/kHaToEV; }
    else if (name=="gridAudit")
    {
        a.kind=PostAction::Kind::GridAudit;
        if (body.contains("reference"))
            throw std::runtime_error("deck: postSCF gridAudit 'reference' (a named raster) is retired with BallOnly/AliasFree (C2 step 3): "
                                     "state 'referenceMargin', the reference rasters' alias-exponent margin over the run's (default "+std::to_string(int(SolidCalculation::kReferenceMargin))+")");
        f.Get("referenceMargin",a.referenceMargin);
    }
    else if (name=="eggBox")
    {
        a.kind=PostAction::Kind::EggBox;
        std::vector<double> v{a.shift.x,a.shift.y,a.shift.z}; f.Get("shift",v);
        if (v.size()!=3) throw std::runtime_error("deck: '"+path+".eggBox.shift' must be [x,y,z] (bohr)");
        a.shift=rvec3_t(v[0],v[1],v[2]);
    }
    else if (name=="xcMeshAudit")
    { a.kind=PostAction::Kind::XCMeshAudit; f.Get("nRadial",a.xcRefNRadial); f.Get("angularDegree",a.xcRefAngularDegree); }
    else if (name=="displace")
    {
        a.kind=PostAction::Kind::Displace;
        const json& pts=f.At("points");
        if (!pts.is_array()) throw std::runtime_error("deck: '"+path+".displace.points' must be a list of {\"atoms\":[i,...], \"delta\":[x,y,z]} (bohr)");
        for (size_t k=0;k<pts.size();++k)
        {
            const std::string pp=path+".displace.points."+std::to_string(k);
            Fields g(pts[k],pp);
            PostAction::Move m;
            g.Get("atoms",m.atoms);
            std::vector<double> v; g.Get("delta",v);
            if (v.size()!=3) throw std::runtime_error("deck: '"+pp+".delta' must be [x,y,z] (bohr)");
            m.delta=rvec3_t(v[0],v[1],v[2]);
            g.Done();
            a.points.push_back(std::move(m));
        }
    }
    else throw std::runtime_error("deck: unknown postSCF action '"+name+"' at "+path+" (legal: estimateHubbardU, hubbardLoop, independentResponse, hubbardLinearResponse, hubbardFiniteDifference, gridAudit, eggBox, xcMeshAudit, displace)");
    f.Done();
    return a;
}

json BasisJson(const RunSpec::Basis& b)
{
    json j={{"data",b.data},{"angular",NameOf(b.angular,kAngularForm)}};
    if (!b.file.empty()) j["file"]=b.file;
    if (!b.inlineBsd.empty()) j["inline"]=b.inlineBsd;
    if (b.autoTrim) j["autoTrim"]=true;
    if (b.family!=RunOptions::Solid::OrbitalBasis::GPW) { j["family"]=NameOf(b.family,kFamily); j["ecut"]=b.ecut; }
    if (!b.trim.empty())
    {
        json t=json::array();
        for (const auto& x : b.trim) t.push_back({{"Z",x.Z},{"l",x.l},{"alpha",x.alpha}});
        j["trim"]=t;
    }
    return j;
}
void BasisFromJson(const json& j, RunSpec::Basis& b)
{
    Fields f(j,"basis");
    f.Get("data",b.data); (void)BasisSet::Gaussian::BasisFileOfName(b.data);   // a name not in the data directory THROWS here, listing the names
    f.GetEnum("angular",b.angular,kAngularForm);
    // DECK-ALIAS: basis.spherical -> basis.angular (2026-11-12)
    if (j.contains("spherical"))   // the earlier spelling (a bool, before 2026-11-12): the same option, so a record still runs
    {
        bool sph=false; f.Get("spherical",sph);
        const auto old=sph ? BasisSet::Gaussian::Angular::Spherical : BasisSet::Gaussian::Angular::Cartesian;
        if (j.contains("angular") && old!=b.angular)
            throw std::runtime_error("deck: basis.angular and the old basis.spherical disagree -- state basis.angular only");
        b.angular=old;
    }
    f.Get("file",b.file); f.Get("inline",b.inlineBsd); f.Get("autoTrim",b.autoTrim);
    if (!b.file.empty() && !b.inlineBsd.empty())
        throw std::runtime_error("deck: basis.file AND basis.inline -- state one (an inline basis IS the file's text)");
    if (j.contains("vet"))
        throw std::runtime_error("deck: 'basis.vet' is retired (2026-11-04): the conditioning report now always runs -- for the vet-stage trim state "
                                 "\"autoTrim\": true beside an orthoTol > 0, and otherwise drop the key");
    f.GetEnum("family",b.family,kFamily); f.Get("ecut",b.ecut);
    if (j.contains("trim"))
    {
        b.trim.clear(); size_t i=0;
        for (const auto& t : f.At("trim"))
        {
            Fields tf(t,"basis.trim."+std::to_string(i++)); RunSpec::Basis::Trim x;
            tf.Get("Z",x.Z); tf.Get("l",x.l); tf.Get("alpha",x.alpha); tf.Done();
            b.trim.push_back(x);
        }
    }
    else f.At("trim");
    f.Done();
}
}

//! A MOLECULE's `basis` (\c RunOptions::Molecule): the same deck key as the solid's, read into the molecule's home -- \c data the
//! basis name (a file stem, checked against the data directory as the solid's is), \c file an override path, \c angular, \c engine.
//! The solid-only keys (inline, autoTrim, trim, family, ecut) are unknown here, so a molecule refuses them by name.
static json MoleculeBasisJson(const RunOptions::Molecule& m)
{
    json j={{"data",m.basis},{"angular",NameOf(m.angular,kAngularForm)},{"engine",NameOf(m.engine,kEngine)}};
    if (!m.basisFile.empty()) j["file"]=m.basisFile;
    return j;
}
static void MoleculeBasisFromJson(const json& j, RunOptions::Molecule& m)
{
    Fields f(j,"basis");
    f.Get("data",m.basis); if (!m.basis.empty()) (void)BasisSet::Gaussian::BasisFileOfName(m.basis);
    f.Get("file",m.basisFile); f.GetEnum("angular",m.angular,kAngularForm); f.GetEnum("engine",m.engine,kEngine);
    f.Done();
}
//! The shared top-level keys a MOLECULE or an ATOM consumes, filled with that kind's defaults (\c MoleculeOptions / \c AtomOptions
//! fill the same values; written here, not by calling them, so writing a record files no policy row).  The keys the kind refuses
//! are not written, so the record reloads to the same run.
static json NonPeriodicJson(const RunOptions& o, RunKind kind)
{
    const bool atom = kind==RunKind::Atom;
    json j={{"spin",NameOf(o.spin,kSpin)},{"model",Hamiltonian::Name(o.model.value_or(Hamiltonian::Model::HF))},{"xalpha",o.xalpha},
            {"seed",NameOf(o.seed,kSeed)},{"accuracy",o.accuracy.eps},{"accelerator",NameOf(o.accelerator,kAcc)},
            {"xc",{{"fit",NameOf(o.vxcFit,kVxc)},{"rho",NameOf(o.xcRho,kRho)},
                   {"mesh",ToJson(o.xcMesh.value_or(atom ? kAtomMeshDefault : kMolecularMeshDefault))}}}};
    if (atom)
    {
        j["ortho"]=NameOf(o.ortho.value_or(qchem::Auto),kOrtho);
        j["atom"]=AtomBlockJson(o.atom);
    }
    else
    {
        json species=json::array();
        for (const auto& s : o.species) species.push_back({s.first,s.second});
        j["coulomb"]={{"fit",NameOf(o.coulomb.fit,kCoulFit)},{"auxBasis",o.coulomb.auxBasis}};
        j.update({{"multiplicity",o.multiplicity},{"pp",o.pp.value_or(std::string())},{"species",species},
                  {"imposeSymmetry",o.imposeSymmetry},{"symmetryTol",o.symmetryTol},{"label",o.label},{"saveStateTo",o.saveStateTo}});
    }
    return j;
}

RunKind KindOf(const RunSpec& r)
{
    if (!r.inlineStructure.is_null()) return RunKind::Solid;   // an inline entry is a materials.json entry: a cell
    if (StructureData::KindOf(r.structure)==StructureData::Kind::Cell) return RunKind::Solid;   // THROWS, listing every name, on a miss
    return StructureData::GetMolecule(r.structure).GetNumAtoms()==1 ? RunKind::Atom : RunKind::Molecule;
}

namespace
{
const char* KindName(RunKind k) { return k==RunKind::Solid ? "a cell" : k==RunKind::Molecule ? "a molecule" : "an atom"; }
}

json ToJson(const RunSpec& r)
{
    const RunKind kind = r.structure.empty() && r.inlineStructure.is_null() ? RunKind::Solid : KindOf(r);   // no structure yet: the solid's layout
    if (kind!=RunKind::Solid)
    {
        json j={{"structure",r.structure},{"schema",kSchemaVersion}};
        if (kind==RunKind::Molecule) j["basis"]=MoleculeBasisJson(r.opts.molecule);
        j.update(NonPeriodicJson(r.opts,kind));
        if (const json a=StatedJson(r.accel); !a.empty()) j["acceleratorOptions"]=a;
        if (!r.state.save.empty() || !r.state.restartFrom.empty())
        {
            json st=json::object();
            if (!r.state.save.empty())        st["save"]=r.state.save;
            if (!r.state.restartFrom.empty()) st["restartFrom"]=r.state.restartFrom;
            j["state"]=st;
        }
        if (r.threads!=ThreadCounts{}) j["threads"]={{"openmp",r.threads.openmp},{"blas",r.threads.blas}};
        if (r.schedule.empty()) j["scf"]=StageScfJson(r.opts.scf,r.opts.recipe);
        else
        {
            json st=json::array();
            for (const auto& s : r.schedule) st.push_back({{"accelerator",NameOf(s.accelerator,kAcc)},{"scf",StageScfJson(s.scf,s.recipe)}});
            j["schedule"]=st;
        }
        return j;
    }
    json j={{"structure", r.inlineStructure.is_null() ? json(r.structure) : r.inlineStructure}};
    j["schema"]=kSchemaVersion;   // a record (and the GUI's resolved deck) is schema 2: the options' structure-neutral keys at the top
    j["kmesh"]={r.kmesh.x,r.kmesh.y,r.kmesh.z}; j["basis"]=BasisJson(r.basis);
    j.update(ToJson(r.opts));    // the run's options: shared keys, accuracy and xc at the top, the periodic ones under `solid`
    j["acceleratorOptions"]=ToJson(ResolveAcceleratorParams(r.accel, SCFAccelerators::AcceleratorParams{}));   // the record: the resolved chain, complete
    if (!r.state.save.empty() || !r.state.restartFrom.empty())
    {
        json st=json::object();
        if (!r.state.save.empty())        st["save"]=r.state.save;
        if (!r.state.restartFrom.empty()) st["restartFrom"]=r.state.restartFrom;
        j["state"]=st;
    }
    if (!r.postSCF.empty()) { json a=json::array(); for (const auto& x : r.postSCF) a.push_back(PostJson(x)); j["postSCF"]=a; }
    if (r.threads!=ThreadCounts{}) j["threads"]={{"openmp",r.threads.openmp},{"blas",r.threads.blas}};   // stated only: absent = the process's
    if (r.schedule.empty()) j["scf"]=StageScfJson(r.opts.scf,r.opts.recipe);    // one stage: the record says it once
    else
    {
        json st=json::array();
        for (const auto& s : r.schedule) st.push_back({{"accelerator",NameOf(s.accelerator,kAcc)},{"scf",StageScfJson(s.scf,s.recipe)}});
        j["schedule"]=st;
    }
    return j;
}
void FromJson(const json& jIn, RunSpec& r)
{
    const json j=RequireSchema2(jIn);   // schema 2 only: a schema-1 deck is refused by name
    Fields f(j,"run");
    if (!j.contains("structure") || !(j.at("structure").is_string() || j.at("structure").is_object()))
        throw std::runtime_error("deck: 'structure' is required: a NAME from materials.json or molecules.json, or the cell itself "
                                 "(an object: lattice, atoms, optional name / pp)");
    if (j.at("structure").is_object()) { r.inlineStructure=f.At("structure"); r.structure.clear(); }
    else                               { f.Get("structure",r.structure); r.inlineStructure=json(); }
    const RunKind kind=KindOf(r);
    if (kind!=RunKind::Solid)   // the periodic-only RunSpec keys, refused by name (the options' own are refused by the kind's view)
    {
        for (const char* k : {"kmesh","postSCF"})
            if (j.contains(k))
                throw std::runtime_error(std::string("deck: '")+k+"' stated for "+KindName(kind)+" ('"+r.structure+"', molecules.json) -- "
                                         +(std::string(k)=="kmesh" ? "it has no Brillouin zone" : "the postSCF actions act on a periodic run")+": drop it");
        if (kind==RunKind::Atom && j.contains("basis"))
            throw std::runtime_error("deck: 'basis' stated for an atom ('"+r.structure+"') -- the radial atom's basis is its 'atom' block "
                                     "(type, basisPool, N / emin / emax, exponents): state it there");
        if (kind==RunKind::Molecule && j.contains("basis")) MoleculeBasisFromJson(f.At("basis"),r.opts.molecule);
        else f.At("basis");
        f.At("kmesh"); f.At("postSCF");
    }
    else if (j.contains("kmesh"))
    {
        const json& k=f.At("kmesh");
        if (!k.is_array() || k.size()!=3) throw std::runtime_error("deck: 'kmesh' must be [nx,ny,nz], got "+k.dump());
        r.kmesh=ivec3_t(k.at(0).get<int>(),k.at(1).get<int>(),k.at(2).get<int>());
    }
    else f.At("kmesh");
    if (kind==RunKind::Solid)
    {
        if (j.contains("basis")) BasisFromJson(f.At("basis"),r.basis); else f.At("basis");
    }
    if (j.contains("state"))
    {
        Fields sf(f.At("state"),"state"); sf.Get("save",r.state.save); sf.Get("restartFrom",r.state.restartFrom); sf.Done();
    }
    else f.At("state");
    if (j.contains("threads"))
    {
        Fields tf(f.At("threads"),"threads"); tf.Get("openmp",r.threads.openmp); tf.Get("blas",r.threads.blas); tf.Done();
        if (r.threads.openmp<0 || r.threads.blas<0)
            throw std::runtime_error("deck: 'threads' counts must be >= 0 (0 = the process's own: QCHEM_OPENMP_THREADS / QCHEM_BLAS_THREADS)");
    }
    else f.At("threads");
    ReadRunOptions(f,j,r.opts);   // the run's options (shared keys, accuracy, xc, the solid block)
    if (j.contains("acceleratorOptions")) FromJson(f.At("acceleratorOptions"),r.accel); else f.At("acceleratorOptions");
    if (j.contains("postSCF"))
    {
        r.postSCF.clear(); size_t i=0;
        for (const auto& a : f.At("postSCF")) r.postSCF.push_back(PostFromJson(a,"postSCF."+std::to_string(i++)));
    }
    else f.At("postSCF");
    if (j.contains("scf") && j.contains("schedule"))
        throw std::runtime_error("deck: give 'scf' (one stage) OR 'schedule' (an annealed recipe), not both -- which would run?");
    if (j.contains("scf")) { FromJson(f.At("scf"),r.opts.scf); r.opts.recipe=RecipeOf(j.at("scf"),"scf"); } else f.At("scf");
    if (j.contains("schedule"))
    {
        r.schedule.clear(); size_t i=0;
        for (const auto& st : f.At("schedule"))
        {
            const std::string path="schedule."+std::to_string(i++);
            Fields sf(st,path); RunSpec::Stage stage;
            sf.GetEnum("accelerator",stage.accelerator,kAcc);
            if (st.contains("scf")) { FromJson(sf.At("scf"),stage.scf); stage.recipe=RecipeOf(st.at("scf"),path+".scf"); } else sf.At("scf");
            sf.Done();
            r.schedule.push_back(stage);
        }
        if (r.schedule.empty()) throw std::runtime_error("deck: 'schedule' is empty");
    }
    else f.At("schedule");
    f.Done();
}
Materials::Material Resolve(RunSpec& spec)
{
    if (!spec.inlineStructure.is_null() && !spec.structure.empty())
        throw std::runtime_error("deck: a structure NAME ("+spec.structure+") AND an inline structure -- state one, not both");
    Materials::Material m;
    if (!spec.inlineStructure.is_null())
    {
        try { m=Materials::Inline(spec.inlineStructure.dump(), spec.opts.pp.value_or("gth-lda-valence")); }
        catch (const std::runtime_error& e) { throw std::runtime_error("deck: structure (inline) -- "+std::string(e.what())); }
    }
    else
    {
        if (StructureData::KindOf(spec.structure)==StructureData::Kind::Molecule)      // KindOf THROWS, listing every known name, on a miss
            throw std::runtime_error("deck: structure '"+spec.structure+"' is a molecule (molecules.json): deck::Prepare builds it, but Run / "
                                     "Vet / View / the basis-shell analysis drive a cell only so far (OpenWork \"Deck: Run / Vet / View for a "
                                     "molecule or an atom\")");
        m=Materials::Get(spec.structure);
    }
    // A NON-SHIPPED BASIS FILE IS INLINED (user ruling 2026-11-10): the elements this material uses, the attribution in the preamble,
    // so the resolved deck -- the revision record -- runs anywhere without the per-user file.
    if (!spec.basis.file.empty() && !ShippedBasisFile(spec.basis.file))
    {
        std::ifstream in(spec.basis.file);
        if (!in) throw std::runtime_error("deck: basis.file '"+spec.basis.file+"' cannot be read");
        const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::vector<std::string> elements;
        m.cell->ForEachSite([&](int z, const rvec3_t&, bool){ elements.push_back(thePeriodicTable().GetSymbol(z)); });
        const auto sub=BasisSet::Gaussian::SubsetBasis(text, elements);
        if (!sub) throw std::runtime_error("deck: basis.file '"+spec.basis.file+"' -- "+sub.Error());
        std::istringstream ls(sub.Value());
        spec.basis.inlineBsd.clear();
        for (std::string l; std::getline(ls, l);) spec.basis.inlineBsd.push_back(l);
        spec.basis.file.clear();
    }
    if (spec.basis.autoTrim && !spec.basis.trim.empty())
        throw std::runtime_error("deck: basis.autoTrim and basis.trim are exclusive (an autoTrim computes the trim; a stated trim is the A/B for it)");
    if (spec.basis.autoTrim && spec.opts.orthoTol<=0.0)
        throw std::runtime_error("deck: basis.autoTrim needs an orthoTol > 0 to trim at -- at 0 it would drop nothing");
    // THE SPHERICAL VIEW ON A CARTESIAN-ONLY BLOCK (OpenWork, 2026-11-14): judged on the FILE the run reads -- the inlined text, a stated
    // basis.file, or the registered basis.data -- by what that file declares (`CARTESIAN-ONLY <elements>`), never by the label in
    // basis.data (which, beside a basis.file, is only a label: a spherical-safe file under the SR label used to be refused, and an SR
    // span under any other label passed).
    if (spec.basis.angular==BasisSet::Gaussian::Angular::Spherical)
    {
        std::set<std::string> cartesianOnly;
        std::string which;
        if (!spec.basis.inlineBsd.empty())
        {
            std::string text; for (const auto& l : spec.basis.inlineBsd) text+=l+"\n";
            std::istringstream is(text);
            cartesianOnly=BasisSet::Gaussian::CartesianOnlyElements(is);
            which="the inlined basis";
        }
        else
        {
            which=spec.basis.file.empty() ? BasisSet::Gaussian::BasisFileOfName(spec.basis.data) : spec.basis.file;
            cartesianOnly=BasisSet::Gaussian::CartesianOnlyElements(which);
        }
        std::set<std::string> hit;
        m.cell->ForEachSite([&](int z, const rvec3_t&, bool)
        {
            std::string sym=thePeriodicTable().GetSymbol(z); for (char& c : sym) c=(char)std::toupper((unsigned char)c);
            if (cartesianOnly.count(sym)) hit.insert(thePeriodicTable().GetSymbol(z));
        });
        if (!hit.empty())
        {
            std::string list; for (const auto& h : hit) list+=(list.empty() ? "" : ", ")+h;
            throw std::runtime_error("deck: basis.angular \"spherical\" with "+which+", whose "+list+" block is CARTESIAN-ONLY (its preamble says so): "
                                     "its s span lives in the Cartesian d contaminants, which the spherical view removes -- use a spherical-safe file "
                                     "(VALENCE_LOWQ_VA, VALENCE_LOWQ_SPH) or basis.angular \"cartesian\"");
        }
    }
    // THE MAGNETIC ORDERING IS STATED (OpenWork §3, 2026-10-26): the facade's own check, run here so a deck is refused BEFORE a revision is claimed.
    try { RequireStatedMagneticOrdering(spec.opts, *m.cell); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: spins -- "+std::string(e.what()).substr(std::string("RunOptions: ").size())); }
    try { RequireConsistentSiteCharges(spec.opts, *m.cell); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: solid.siteCharges -- "+std::string(e.what()).substr(std::string("RunOptions: ").size())); }
    // AUTO PICKS THE KIND, NOT THE RECIPE (OpenWork §3, 2026-10-26): the mesh policy's own refusal, run here so --check catches it too.
    try { (void)qcMesh::ResolveXCMesh(spec.opts.xcMesh.value_or(qcMesh::XCMeshDefault())); }
    catch (const std::invalid_argument& e) { const std::string w=e.what(); throw std::runtime_error("deck: xc.mesh"+(w.rfind("xcMesh",0)==0 ? w.substr(6) : ": "+w)); }

    // PRE-FLIGHT of the postSCF list (each message names the action and what to change): the old probes found these out AFTER the SCF.
    {
        size_t i=0;
        bool displaced=false;
        for (const auto& a : spec.postSCF)
        {
            const std::string who="deck: postSCF."+std::to_string(i++)+" ("+PostName(a.kind)+"): ";
            if (a.kind==PostAction::Kind::GridAudit)   // the one action that is not about +U: it needs a raster, and a different one
            {
                if (spec.basis.family!=RunOptions::Solid::OrbitalBasis::GPW) throw std::runtime_error(who+"a plane-wave run has no raster to audit");
                if (!(a.referenceMargin>0.0)) throw std::runtime_error(who+"referenceMargin must be > 0 (the reference is FINER than the run)");
                continue;
            }
            if (a.kind==PostAction::Kind::XCMeshAudit) // the run's converged D on a reference Becke mesh: a Gaussian run with a positive mesh
            {
                if (spec.basis.family!=RunOptions::Solid::OrbitalBasis::GPW) throw std::runtime_error(who+"a plane-wave run samples XC on its FFT grid, not a Becke mesh");
                if (a.xcRefNRadial<=0 || a.xcRefAngularDegree<=0) throw std::runtime_error(who+"nRadial and angularDegree must be > 0");
                // displace runs in its own session, so this run's D stays on the deck's positions: after one, the audit would read as the displaced point's.
                if (displaced) throw std::runtime_error(who+"it audits the deck's own positions, not a displaced point -- list it before displace, or state the displaced structure in the deck");
                continue;
            }
            if (a.kind==PostAction::Kind::EggBox)      // a rigid translation of a free, U = 0 Gaussian run at a nonzero shift
            {
                if (spec.basis.family!=RunOptions::Solid::OrbitalBasis::GPW) throw std::runtime_error(who+"a plane-wave basis is translation-invariant by construction");
                if (spec.opts.imposeSymmetry) throw std::runtime_error(who+"an imposed group's SALC phases follow the origin -- gate the free twin (imposeSymmetry false)");
                if (!spec.opts.solid.hubbard.empty()) throw std::runtime_error(who+"the +U projectors would not follow the shift -- gate the U = 0 twin");
                if (norm(a.shift)==0.0) throw std::runtime_error(who+"a zero shift measures nothing");
                continue;
            }
            if (a.kind==PostAction::Kind::Displace)    // geometry points through ONE SolidSession (its own setup, sized on max|dH_ij|)
            {
                if (a.points.empty()) throw std::runtime_error(who+"no points -- list {\"atoms\":[i,...], \"delta\":[x,y,z]} (bohr)");
                displaced=true;
                // The session's own refusal, run here so a deck fails before its SCF: a free mesh made invariant under the group detected at
                // EACH geometry switches along a path that breaks it (OpenWork, closed 2026-10-30).
                if (!spec.opts.imposeSymmetry && spec.basis.family==RunOptions::Solid::OrbitalBasis::GPW && spec.opts.solid.xcMeshSymmetric)
                    throw std::runtime_error(who+"a free Gaussian run's XC mesh must be the plain one along a path -- set solid.xcMeshSymmetric false");
                for (size_t k=0;k<a.points.size();++k)
                {
                    const auto& mv=a.points[k];
                    if (norm(mv.delta)==0.0) throw std::runtime_error(who+"points."+std::to_string(k)+": a zero delta is the reference itself");
                    for (size_t at : mv.atoms)
                        if (at>=m.cell->GetNumAtoms())
                            throw std::runtime_error(who+"points."+std::to_string(k)+": atom "+std::to_string(at)+" is past the structure's "
                                                     +std::to_string(m.cell->GetNumAtoms())+" atoms");
                }
                continue;
            }
            const bool needsU = true;   // every action works on the Hubbard manifolds (list them at U=0 to probe a run without +U)
            if (needsU && spec.opts.solid.hubbard.empty())
                throw std::runtime_error(who+"the run carries no Hubbard manifold -- list the channels in solid.hubbard (at U_eV=0 to probe without +U)");
            const bool lrLoop = a.kind==PostAction::Kind::HubbardLoop && a.estimator==SolidCalculation::HubbardLoop::Estimator::LinearResponse;
            const bool isResponse = a.kind==PostAction::Kind::IndependentResponse || a.kind==PostAction::Kind::HubbardLinearResponse || lrLoop;
            // NOT an error: whether an imposed group REDUCES the k-mesh depends on the group and the mesh (the NiO AFM-II 2x2x2 Shubnikov group maps 8 k-points onto 8
            // irreducible ones, and its chi0 ran).  A blanket refusal would reject decks the library accepts; the response itself THROWS at run time if the mesh was
            // really reduced (D5).  Vet() says so in a note.
            // A smeared (metal) response at q != 0 needs ONE Fermi level across the k-mesh (Reference::Gap refuses it at run time --
            // after the SCF; SrVO3 2026-10-05 ran two hours to find that out).  The FINAL stage's kT is the response's occupancy rule.
            if (isResponse && spec.Stages().back().scf.SmearingkT>0.0 && spec.opts.solid.perKFermiLevels)
            {
                const ivec3_t q = a.nq;
                if (q.x*q.y*q.z>1)
                    throw std::runtime_error(who+"a SMEARED final stage (kT "+std::to_string(spec.Stages().back().scf.SmearingkT)+" Ha) with per-k "
                                             "Fermi levels: a metal's response at q != 0 needs ONE Fermi level -- drop solid.perKFermiLevels");
            }
            // littleGroupK MEASURES the ground state's group on the reference's orbitals to 1e-6 (KAxisPlan §5 ruling 1): a reference
            // stopped at Δρ/N 1e-6 sits AT that tolerance and keeps a random subset of ops (TiO2 k221 2026-10-28: 4..16 of 16 at 1e-6;
            // 16 of 16, worst 9e-9, at 1e-9).  Not floored by ρ_lost/N: Δρ is a difference on ONE quadrature, the charge error cancels.
            if (a.kind==PostAction::Kind::HubbardLinearResponse && a.littleGroupK && LooserThan(spec.Stages().back().scf.MinΔρPerElectron, kResponseReferenceMinΔρ))
                throw std::runtime_error(who+"littleGroupK measures the reference's symmetry to 1e-6, so the final SCF stage must converge to "
                                         "scf.minDeltaRhoPerElectron <= 1e-9 (it is "+GateJson(spec.Stages().back().scf.MinΔρPerElectron).dump()+")");
            // (No complex demand: the setup decides the ansatz ONCE from k and q -- PlannedResponseQMeshes, below -- and a real
            //  ground state hands a response its state through the two-run split, user ruling 2026-11-18.)
            if (isResponse)             // the q points loop INSIDE the response, on the k-mesh's own differences
            {
                const ivec3_t& q=a.nq;
                if (q.x<1 || q.y<1 || q.z<1) throw std::runtime_error(who+"every nq division must be >= 1");
                if (spec.kmesh.x%q.x || spec.kmesh.y%q.y || spec.kmesh.z%q.z)
                    throw std::runtime_error(who+"the q-mesh "+QMeshJson(q).dump()+" is incommensurate with kmesh ["+std::to_string(spec.kmesh.x)+","
                                             +std::to_string(spec.kmesh.y)+","+std::to_string(spec.kmesh.z)+"] (each k division must be a multiple of nq's)");
            }
            if (a.kind==PostAction::Kind::HubbardFiniteDifference && a.alpha==0.0)
                throw std::runtime_error(who+"alpha_eV must be nonzero (the +-alpha step)");
            if (a.kind==PostAction::Kind::HubbardLoop && (a.maxOuter<1 || a.tolU_eV<=0.0))
                throw std::runtime_error(who+"needs maxOuter >= 1 and tolU_eV > 0");
            if ((a.kind==PostAction::Kind::HubbardLinearResponse || lrLoop) && (a.maxIter<1 || a.restart<1 || a.tol<=0.0))
                throw std::runtime_error(who+"needs maxIter >= 1, restart >= 1 and tol > 0");
            for (size_t j : a.perturb)
                if (j>=spec.opts.solid.hubbard.size())
                    throw std::runtime_error(who+"perturb index "+std::to_string(j)+" is past the "+std::to_string(spec.opts.solid.hubbard.size())+" manifolds in solid.hubbard");
        }
    }
    // The material NAMES its pseudopotential set; a deck that left `pp` at the facade default takes the material's.  `species` (the
    // per-element overrides) is filled with the resolved vocabulary so the record states every valence actually used.
    if (!spec.opts.pp) spec.opts.pp=m.pp;   // unstated: the material's set
    if (spec.opts.species.empty()) spec.opts.species=m.species;
    // The electron count of the species the RUN uses: a stated override (Na q9) changes it.  Until 2026-10-31 this read the
    // material's own species, so a deck stating Na q9 ran NaF with 8 electrons instead of 16 -- a cell charged +8, silently.
    if (spec.opts.Nelec==0)
    {
        Materials::Material used=m;
        for (auto& [el,q] : used.species)
            for (const auto& [el2,q2] : spec.opts.species) if (el2==el) q=q2;
        spec.opts.Nelec=used.Nelec()-spec.opts.charge;   // a STATED net charge (RunOptions::charge)
    }
    // THE SPIN BOOKKEEPING on the resolved count (the setup's own decision 2, run here so --check refuses an impossible request:
    // multiplicity 1 on an odd count once ran a silent doublet).
    try { (void)SpinChannelCounts(spec.opts, spec.opts.Nelec); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: multiplicity -- "+std::string(e.what()).substr(std::string("RunOptions: ").size())); }
    // THE MIXER'S OWN ADMISSION (the iterator's ResolveMixing), run here so --check refuses what the run would: Pulay with no pulayDepth,
    // Kerker with no kerkerG0, a density history beside a STATED Fock extrapolator (DIIS, Ladder).  A deck run is periodic (gSpace).  An
    // Auto accelerator is judged as the cold-stage policy's Null: on a RESTORED stage the run may pick DIIS / Ladder (MakeAccelerator), and
    // only the run can say so -- the one case --check leaves to it.
    {
        using A=SCFAccelerators::Type;
        const auto stages=spec.Stages();
        for (size_t i=0; i<stages.size(); ++i)
        {
            const bool fock = stages[i].accelerator==A::DIIS || stages[i].accelerator==A::Ladder;
            const std::string where = spec.schedule.empty() ? std::string("scf") : "schedule."+std::to_string(i)+".scf";
            // THE RECIPE'S ADMISSION (qchem.SCFRecipe): MagneticInsulator / Metal refuse kT = 0.  Then the mixer as the recipe fills it.
            try { AdmitRecipe(stages[i].recipe, stages[i].scf.SmearingkT, /*molecule*/false); }
            catch (const std::invalid_argument& e) { throw std::runtime_error("deck: "+where+".recipe -- "+e.what()); }
            const SCFStageKind kind = i==0 ? SCFStageKind::Cold : SCFStageKind::Carried;
            try { (void)qchem::ResolveMixing(FillRecipeKnobs(stages[i].recipe, kind, stages[i].scf).params, /*gSpace*/true, fock); }
            catch (const std::invalid_argument& e)
            {
                throw std::runtime_error("deck: "+where+".mixer -- "+e.what());
            }
        }
    }
    return m;
}

//=== --set =======================================================================================
void ApplySet(json& deck, std::string_view assignment)
{
    const size_t eq=assignment.find('=');
    if (eq==std::string_view::npos || eq==0) throw std::runtime_error("deck: --set needs path=value, got '"+std::string(assignment)+"'");
    const std::string path(assignment.substr(0,eq)), valueText(assignment.substr(eq+1));
    json value = json::parse(valueText, nullptr, /*exceptions*/false);
    if (value.is_discarded()) value=valueText;                       // a bare word is a string
    json* at=&deck;
    for (size_t b=0; b<=path.size();)
    {
        const size_t e=std::min(path.find('.',b), path.size());
        const std::string key=path.substr(b,e-b);
        if (key.empty()) throw std::runtime_error("deck: --set path '"+path+"' has an empty component");
        const bool last=e>=path.size();
        if (at->is_array())
        {
            const size_t idx=std::stoul(key);
            if (idx>=at->size()) throw std::runtime_error("deck: --set '"+path+"': index "+key+" is past the end of an array of "+std::to_string(at->size()));
            at=&(*at)[idx];
        }
        else
        {
            if (!at->is_object()) *at=json::object();
            at=&(*at)[key];
        }
        if (last) { *at=value; return; }
        b=e+1;
    }
}

//=== revisions ===================================================================================
std::filesystem::path ClaimRevision(const std::filesystem::path& dir, const std::string& name)
{
    std::filesystem::create_directories(dir);
    for (int n=1; n<=99999; ++n)
    {
        char num[16]; std::snprintf(num,sizeof num,"r%03d",n);
        const std::filesystem::path p=dir/(name+"."+num+".json");
        const int fd=::open(p.c_str(), O_WRONLY|O_CREAT|O_EXCL, 0644);   // atomic: exactly one claimant wins a number
        if (fd>=0) { ::close(fd); return p; }
    }
    throw std::runtime_error("deck: no free revision number for '"+name+"' in "+dir.string());
}

namespace
{
std::string Fnv1a(const std::filesystem::path& p)   // a checksum of the input deck, for the provenance block (identity, not security)
{
    std::ifstream in(p,std::ios::binary); unsigned long long h=1469598103934665603ULL; char c;
    while (in.get(c)) { h^=(unsigned char)c; h*=1099511628211ULL; }
    char b[32]; std::snprintf(b,sizeof b,"%016llx",h); return b;
}
}

void WriteRevision(const std::filesystem::path& path, const json& resolved, const Provenance& prov)
{
    json out;
    out["deck"]={{"schema",kSchemaVersion},{"codeVersion",prov.codeVersion}};
    json pv={{"commandLine",prov.commandLine},{"overrides",prov.overrides},{"ignoredEnvironment",prov.ignoredEnvironment},
                   {"activeEnvironment",prov.activeEnvironment}};
    if (!prov.policyResolved.empty()) pv["policyResolved"]=prov.policyResolved;
    if (!prov.restartedFrom.empty())  pv["restartedFrom"]=prov.restartedFrom;
    if (!prov.basisFile.empty())      { pv["basisFile"]=prov.basisFile; pv["basisChecksum"]=prov.basisChecksum; }
    if (!prov.inputDeck.empty())
    {
        pv["inputDeck"]=prov.inputDeck.string();
        if (std::filesystem::exists(prov.inputDeck))
        {
            pv["inputDeckChecksum"]=Fnv1a(prov.inputDeck);
            std::ifstream in(prov.inputDeck); const json src=json::parse(in,nullptr,false);
            if (src.is_object() && src.contains("deck") && src.contains("run")) pv["parent"]=prov.inputDeck.filename().string();   // started from a revision: the lineage
        }
    }
    out["provenance"]=pv;
    out["run"]=resolved;
    std::ofstream os(path,std::ios::trunc);
    if (!os) throw std::runtime_error("deck: cannot write "+path.string());
    os<<out.dump(2)<<"\n";
}

void ResolveBasisFile(json& payload, const std::filesystem::path& deckDir)
{
    if (!payload.contains("basis") || !payload["basis"].contains("file")) return;
    const std::filesystem::path f=payload["basis"]["file"].get<std::string>();
    if (f.is_relative() && std::filesystem::exists(deckDir/f)) payload["basis"]["file"]=(deckDir/f).string();
}

json LoadDeck(const std::filesystem::path& path, const std::string& currentCodeVersion)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("deck: cannot read "+path.string());
    json j;
    try { j=json::parse(in); } catch (const json::exception& e) { throw std::runtime_error("deck: "+path.string()+": "+e.what()); }
    if (!j.contains("run")) { ResolveBasisFile(j,path.parent_path()); return j; }   // a hand-written deck: the payload IS the file
    const json& h=j.value("deck",json::object());
    if (h.value("schema",0) > kSchemaVersion)
        throw std::runtime_error("deck: "+path.string()+" has schema "+std::to_string(h.value("schema",0))+", newer than this code's "+std::to_string(kSchemaVersion));
    if (h.value("codeVersion",std::string()) != currentCodeVersion)
        std::fprintf(stderr,"[deck] WARNING: %s was written by code version '%s', this is '%s' -- the run may not reproduce\n",
                     path.c_str(), h.value("codeVersion",std::string()).c_str(), currentCodeVersion.c_str());
    json payload=j.at("run");
    ResolveBasisFile(payload,path.parent_path());
    return payload;
}


static std::string Refusal(const SetupRefusal& r);   // the one wording of a setup's refusal (defined beside Prepare)

VetReport Vet(RunSpec spec, const std::filesystem::path& outDir)
{
    namespace fs=std::filesystem;
    // A MOLECULE or an ATOM deck (Group face step 3b; controller 2026-11-21: the early pin-30 refusal end to end): its own setup's Admit --
    // every policy decision, the Vxc pre-flight among them, and nothing built (no basis, no mesh) -- refused with the setup's own reason.
    if (const RunKind kind=KindOf(spec); kind!=RunKind::Solid)
    {
        const auto& st=StructureData::GetMolecule(spec.structure);
        if (kind==RunKind::Molecule)
        {
            if (auto s=MoleculeSetup::Admit(st, spec.opts, spec.accel); !s) throw std::runtime_error(Refusal(s.Error()));
        }
        else if (auto s=AtomSetup::Admit(st, spec.opts, spec.accel); !s) throw std::runtime_error(Refusal(s.Error()));
        VetReport v;
        v.summary=spec.StructureName()+": "+(kind==RunKind::Molecule ? "molecule" : "atom")+", its setup admitted (no basis or mesh built)";
        return v;
    }
    const Materials::Material mat=Resolve(spec);
    // A GGA: Run's setup admits it on the free Becke delta route only (GGARoute, GGA stage 3), so Vet asks the same rule of what the deck
    // STATES.  An unstated (Auto) mesh is what the setup then CHOOSES for a GGA, Becke.  (OpenWork G4: Vet should run the setup itself.)
    if (spec.opts.model)
    {
        qcMesh::MeshParams mesh=spec.opts.xcMesh.value_or(qcMesh::XCMeshDefault());
        if (mesh.cellKind==qcMesh::UnitCellKind::Auto) mesh.cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Becke);   // what the setup chooses (above)
        RouteSolidGradient(spec.opts, mesh, spec.opts.imposeSymmetry);
    }
    // Pin 30, the stated Vxc fit and XC mesh on this cell (Group face step 3b, the user's condition on ruling 2026-11-21: Vet keeps refusing
    // a bad deck before Run).  The setup's own cheap pre-flight, asked of the cell's group alone: nothing is built.  The STATED mesh is
    // asked: on a cell its Auto resolves to Becke or uniform, and the rule gives both the same answer.
    try
    {
        RunOptions o=spec.opts;
        RequireVxcAdmitted(*mat.cell, o.vxcFit, SolidXCMesh(o), "cell (deck)");
    }
    catch (const policy::Refused& r) {throw std::runtime_error(Refusal(SetupRefusal{r.Name(), r.what()}));}   // the setups' one wording
    VetReport v;
    if (!spec.state.restartFrom.empty())
    {
        const std::string r=spec.state.restartFrom;
        const bool isPath = r.find('/')!=std::string::npos || (r.size()>3 && r.substr(r.size()-3)==".h5");
        const std::string path = isPath ? r : (outDir/"states"/(r+".h5")).string();
        if (!fs::exists(path)) v.notes.push_back("state.restartFrom: no saved state at "+path+" (a run would refuse)");
    }
    Lattice_3D lat(*mat.cell, spec.kmesh);
    GaussianBasisRecipe recipe=RecipeOf(spec.basis); recipe.autoTrim=false;   // --check builds the stated basis; the vet loop is a run's
    std::shared_ptr<const BasisSet::Real_BS> b=recipe.Make(*mat.cell, Pseudopotential::PseudopotentialSet::Get(spec.opts.pp.value_or("gth-lda-valence")).With(spec.opts.species));
    size_t nf=0;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*b).Iterate<BasisSet::Real_OIBS>()) nf+=ibs->GetNumFunctions();
    {
        bool resp=false; for (const auto& a : spec.postSCF) resp = resp || a.kind==RunSpec::PostAction::Kind::IndependentResponse || a.kind==RunSpec::PostAction::Kind::HubbardLinearResponse;
        if (resp && spec.opts.imposeSymmetry && spec.kmesh.x*spec.kmesh.y*spec.kmesh.z>1)
            v.notes.push_back("a response with imposeSymmetry on a k-mesh: if the imposed group REDUCES the mesh the response throws at run time (D5); the NiO AFM-II 2x2x2 Shubnikov group does not reduce it");
        if (resp && LooserThan(spec.Stages().back().scf.MinΔρPerElectron, kResponseReferenceMinΔρ))
            v.notes.push_back("a response reference should converge to scf.minDeltaRhoPerElectron <= 1e-9 (a few more iterations; ρ_lost/N does not floor it): "
                              "at 1e-6 its orbitals carry ~1e-6 symmetry and eigenpair residuals");
    }
    if (spec.basis.autoTrim) v.notes.push_back("basis.autoTrim: the vet-stage trim loop is NOT run by --check (function count is the untrimmed one)");
    std::ostringstream os;
    os<<spec.StructureName()<<(spec.inlineStructure.is_null() ? "" : " (inline)")<<": "<<mat.cell->GetNumAtoms()<<" atoms, Nelec "<<spec.opts.Nelec<<", multiplicity "<<spec.opts.multiplicity
      <<", k-mesh "<<spec.kmesh.x<<"x"<<spec.kmesh.y<<"x"<<spec.kmesh.z<<", basis "<<(!spec.basis.inlineBsd.empty() ? spec.basis.data+" (inline)" : spec.basis.file.empty()?spec.basis.data:spec.basis.file)<<(spec.basis.angular==BasisSet::Gaussian::Angular::Spherical?" (spherical)":"")
      <<" "<<nf<<" functions, "<<spec.Stages().size()<<" SCF stage(s), "<<spec.opts.solid.hubbard.size()<<" Hubbard manifold(s), "<<spec.postSCF.size()<<" postSCF action(s)"
      <<(spec.opts.imposeSymmetry?", symmetry IMPOSED":", symmetry FREE");
    v.summary=os.str();
    return v;
}

// postSCF displace: the deck's own setup (re-resolved: an UNSIZED one, so the session's reference pilot sizes it on max|dH_ij|), the
// reference at the deck's positions, then each point -- all through ONE SolidSession, so every point shares the reference's frozen
// basis recipe and rasters.  The main run is not reused: it was sized as an Energy.  No saved state (the main run's file is the deck's).
//! WHICH RUN of a deck a setup is admitted for (the two-run split, user ruling 2026-11-18).  The ansatz is decided ONCE per run, from k
//! and q: \c Whole and \c Response state the deck's planned linear responses (every block complex, \c SolidSetup::realTRIMBlocks);
//! \c GroundState states none, so its TRIM blocks run real.
enum class RunPart {Whole, GroundState, Response};
static std::shared_ptr<const SolidSetup> SetupOf(const RunSpec&, const Materials::Material&, const std::vector<SCFStage>&,
                                                 RunPart part=RunPart::Whole);
// Build's body with the setup's refusal as a value (Run reports it, Build throws it); defined beside Build.
static Outcome<std::unique_ptr<SolidCalculation>, SetupRefusal> PrepareOf(const RunSpec& specIn, SolidCalculation::Observer onIteration,
                                                                          RunPart part=RunPart::Whole);
// Run #2 of the split: the deck's RESPONSE setup, the saved ground state at \a statePath adopted on it (SolidCalculation::Adopt).
static std::unique_ptr<SolidCalculation> AdoptOf(const RunSpec& specIn, const std::string& statePath);

//! The q-meshes of the deck's planned LINEAR responses (a hubbardLinearResponse, a linear-response U loop, and a hubbardFiniteDifference
//! as q = 0: its +-alpha SCFs measure occupations through the complex-only Hubbard probe), once each: what the setup decides the ansatz
//! from.  (An independentResponse is not one: chi0 runs on real TRIM blocks.)
static std::vector<ivec3_t> PlannedResponseQMeshes(const RunSpec& spec)
{
    std::vector<ivec3_t> out;
    for (const auto& a : spec.postSCF)
    {
        const bool lrLoop = a.kind==PostAction::Kind::HubbardLoop && a.estimator==SolidCalculation::HubbardLoop::Estimator::LinearResponse;
        const bool fd     = a.kind==PostAction::Kind::HubbardFiniteDifference;
        if (a.kind!=PostAction::Kind::HubbardLinearResponse && !lrLoop && !fd) continue;
        const ivec3_t q = fd ? ivec3_t(1,1,1) : a.nq;
        bool known=false;
        for (const ivec3_t& p : out) known = known || (p.x==q.x && p.y==q.y && p.z==q.z);
        if (!known) out.push_back(q);
    }
    return out;
}
//! \brief Does the deck run as TWO runs (user ruling 2026-11-18)?  Run #1, the SCF and every action before the first linear response,
//! on a setup with no response planned (real TRIM blocks: the SCF ~1.3x faster and lighter, TiO₂ k222 1619 vs 2139 s); run #2, the
//! responses, on a complex setup that ADOPTS run #1's converged state (no SCF).  Only when it can work: a linear response is planned,
//! no U loop re-runs a complex SCF (a linear-response loop), no finite-difference chi runs its +-alpha SCFs (on a complex setup, from
//! wherever it sits), and nothing after the first response needs an SCF'd run (only responses follow it).  Otherwise the deck is one
//! run, complex from the start when a response is planned.
static bool SplitsAtResponse(const RunSpec& spec)
{
    if (PlannedResponseQMeshes(spec).empty()) return false;
    bool inResponses=false;
    for (const auto& a : spec.postSCF)
    {
        if (a.kind==PostAction::Kind::HubbardLoop && a.estimator==SolidCalculation::HubbardLoop::Estimator::LinearResponse) return false;
        if (a.kind==PostAction::Kind::HubbardFiniteDifference) return false;
        const bool response = a.kind==PostAction::Kind::HubbardLinearResponse || a.kind==PostAction::Kind::IndependentResponse;
        if (a.kind==PostAction::Kind::HubbardLinearResponse) inResponses=true;
        if (inResponses && !response) return false;
    }
    return true;
}
struct DisplaceResult
{
    struct Point { std::optional<double> energy; size_t iterations=0; std::string failure; };
    std::optional<double> E0; size_t it0=0; std::string failure0;
    std::vector<Point> points;
    bool ok=false;
};
static DisplaceResult Displace(const RunSpec& specIn, const std::vector<RunSpec::PostAction::Move>& moves)
{
    RunSpec spec=specIn;
    const Materials::Material mat=Resolve(spec);
    spec.opts.solid.orbitalBasis=spec.basis.family; spec.opts.solid.planeWaveEcut=spec.basis.ecut;
    spec.opts.saveStateTo.clear();
    std::vector<SCFStage> stages;
    for (const auto& st : spec.Stages()) stages.push_back({st.scf,st.accelerator,st.recipe});   // the whole recipe: the session seeds, it never restarts
    SolidSession S(SetupOf(spec, mat, stages), stages);
    const PESGeometry ref=S.Reference();
    DisplaceResult out;
    std::cout<<"[displace] reference: the deck's positions ("<<ref.R.size()<<" atoms), one SolidSession for "<<moves.size()<<" point(s)"<<std::endl;
    // Each point in a row of its own: the point writes its timing breakdown (and anything else it reports) one level in.
    auto p0=[&]{ const qchem::report::Row row("displacePoints"); return S.ConvergeAt(ref.cell, ref.R); }();
    if (!p0) { out.failure0=p0.Error().details; return out; }
    out.E0=p0->energy; out.it0=p0->iterations;
    bool ok=true;
    for (size_t k=0;k<moves.size();++k)
    {
        std::vector<rvec3_t> R=ref.R;
        if (moves[k].atoms.empty()) for (auto& r : R) r+=moves[k].delta;
        else                        for (size_t a : moves[k].atoms) R[a]+=moves[k].delta;
        auto p=[&]{ const qchem::report::Row row("displacePoints"); return S.ConvergeAt(ref.cell, R); }();
        DisplaceResult::Point q;
        if (p) { q.energy=p->energy; q.iterations=p->iterations; }
        else   { q.failure=p.Error().details; ok=false; }
        std::ostringstream os;
        if (p) os<<"E "<<std::setprecision(12)<<p->energy<<" Ha, dE "<<std::scientific<<std::setprecision(6)<<p->energy-*out.E0<<" Ha, "<<p->iterations<<" iterations";
        else   os<<"FAILED: "<<q.failure;
        std::cout<<"[displace] point "<<k<<": "<<os.str()<<std::endl;
        out.points.push_back(std::move(q));
    }
    out.ok=ok;
    return out;
}

RunOutcome Run(RunSpec spec, Provenance prov, const std::filesystem::path& outDir)
{
    namespace fs=std::filesystem;
    struct Close { ~Close() { qchem::report::ClearConsole(); qchem::report::End(); } };
    // A run is OPEN before anything can refuse, so a refusal files its record (its policy rows) too.  The checks before the
    // revision is claimed (Resolve, the restart) run in "<structure> resolve"; the build and SCF in the revision's own run.
    const qchem::ScopedWorkerThreads threads(spec.threads);   // the run's thread counts, for the whole run (absent = the process's)
    std::optional<Close> resolving;
    qchem::report::Begin(spec.StructureName()+" resolve");
    qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
    resolving.emplace();
    (void)Resolve(spec);   // the record states the derived Nelec / species; Build resolves its own copy

    // RESTART: resolve (and check) the saved state BEFORE claiming a revision, so a typo leaves no half-record.  A bare stem names a revision of this
    // outDir (-> its states/<stem>.h5, and it becomes this run's lineage); anything with a '/' or an .h5 suffix is a path.
    std::string restartPath;
    if (!spec.state.restartFrom.empty())
    {
        const std::string r=spec.state.restartFrom;
        const bool isPath = r.find('/')!=std::string::npos || (r.size()>3 && r.substr(r.size()-3)==".h5");
        restartPath = isPath ? r : (outDir/"states"/(r+".h5")).string();
        if (!fs::exists(restartPath)) throw std::runtime_error("deck: state.restartFrom '"+r+"': no saved state at "+restartPath);
        if (!isPath) prov.restartedFrom=r+".json";
    }

    resolving.reset();
    RunOutcome out;
    out.revision=ClaimRevision(outDir, spec.StructureName());
    prov.policyResolved=RunPolicy(spec.opts.solid.policy).Banner();   // the stated policy is in the deck; what it RESOLVED to is the record
    prov.basisFile=!spec.basis.inlineBsd.empty() ? MaterializedBasis(spec.basis.inlineBsd)
                  : spec.basis.file.empty() ? BasisSet::Gaussian::BasisFileOfName(spec.basis.data) : spec.basis.file;
    prov.basisChecksum=BasisSet::Gaussian::BasisFileChecksum(prov.basisFile);
    WriteRevision(out.revision, ToJson(spec), prov);     // the record exists BEFORE the SCF: a crashed run still leaves it
    const std::string stem=out.revision.stem().string(); // <structure>.rNNN -- the run's name in its own output

    if (spec.opts.label=="gpw") spec.opts.label=stem;
    if (!spec.state.save.empty())
    {
        if (spec.state.save=="auto") { fs::create_directories(outDir/"states"); spec.state.save=(outDir/"states"/(stem+".h5")).string(); }
        // A STATED path's directory too: a missing one made the save WARN after the whole SCF and the run read CONVERGED with no state
        // written (OpenWork §3, TiO₂ 2026-11-18: a 26-min run lost).
        else if (const fs::path dir=fs::path(spec.state.save).parent_path(); !dir.empty()) fs::create_directories(dir);
    }
    if (!restartPath.empty()) spec.state.restartFrom=restartPath;   // Build takes a path: the stem was this outDir's
    const std::vector<SCFStage> stages=ScheduleToRun(spec);

    qchem::report::Begin(stem);
    qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
    const Close close;
    // A DRIVER timing window for the whole deck run: its calculations stay quiet (the deck emits ONE `timing` section at its end, below,
    // what scripts/bench greps); a displaced point or a geometry step still reports its own.
    qchem::report::TimingWindow deckTiming(stem);
    // THE TWO-RUN SPLIT (SplitsAtResponse): run #1 is the ground state on a setup with no response planned (real TRIM blocks); the
    // responses run on run #2, which adopts run #1's state (below, at the first linear response).
    const bool split=SplitsAtResponse(spec);
    auto prepared=PrepareOf(spec, nullptr, split ? RunPart::GroundState : RunPart::Whole);   // ONE path with the GUI's (Build): built, then converged here
    if (!prepared)
    {   // REFUSED by the setup: the outcome says so and the record (claimed above) files it -- no SCF ran, so no postSCF action did either
        out.summary="REFUSED (policy "+prepared.Error().policy+"): "+prepared.Error().reason;
        std::ifstream in(out.revision); json rec=json::parse(in); in.close();
        rec["results"]={{"converged",false},{"refused",{{"policy",prepared.Error().policy},{"reason",prepared.Error().reason}}}};
        std::ofstream os(out.revision,std::ios::trunc); os<<rec.dump(2)<<"\n";
        return out;
    }
    std::unique_ptr<SolidCalculation> calc=prepared.TakeValue();   // not const: the split replaces run #1 by run #2
    (void)calc->Converge(stages);
    if (auto R=calc->Result())
    {
        out.converged=true; out.energy=calc->LastIterateTerms().GetTotalEnergy();
        out.summary="CONVERGED  "+calc->Diagnostics().Summary();
    }
    else out.summary=(R.Error().why==SCFFailure::Why::ResourceRefused ? "REFUSED ON RESOURCES: "
                    : R.Error().why==SCFFailure::Why::StageRefused    ? "REFUSED: "
                    : R.Error().why==SCFFailure::Why::StateNotSaved   ? "converged, STATE NOT SAVED: " : "NOT converged: ")+R.Error().details;
    // A run REFUSED ON RESOURCES (the audit RSS policy) did not reach the state its postSCF actions read: before the first SCF there is no D
    // at all, and mid-pilot the rasters are not settled.  Each action is recorded as skipped, and the record carries the refusal's numbers.
    std::optional<ResourceRefusal> refused;
    if (const auto R=calc->Result(); !R && R.Error().resource) refused=R.Error().resource;

    // postSCF: each action on the converged calculation, in order, each reporting itself; the summaries go in the revision's `results`.  A failed
    // action (a response gate, an unconverged FD step) is a RESULT, not an exception: the remaining actions still run.
    json results=json::array();
    bool handedOver=false;   // the split's run #2 has replaced run #1
    for (const auto& a : spec.postSCF)
    {
        RunOutcome::PostResult pr; pr.action=PostName(a.kind);
        json rj={{"action",pr.action}};
        if (refused)
        {
            pr.ok=false; pr.summary="SKIPPED: the run was refused on resources ("+refused->when+")";
            std::cout<<"["<<stem<<"] postSCF: "<<pr.action<<"  "<<pr.summary<<std::endl;
            rj["ok"]=pr.ok; rj["summary"]=pr.summary; results.push_back(rj);
            out.postSCF.push_back(pr);
            continue;
        }
        std::cout<<"["<<stem<<"] postSCF: "<<pr.action<<(out.converged ? "" : "  (the SCF did NOT converge -- a diagnostic only)")<<std::endl;
        using K=PostAction::Kind;
        // THE HAND-OVER (the two-run split, user ruling 2026-11-18): at the first linear response, run #1's converged state goes to
        // file, run #1 is DESTROYED (its real tables with it), and run #2 -- the response setup, every block complex -- adopts the state
        // with one Fock check and no SCF (SolidCalculation::Adopt).  Every action from here on is a response (SplitsAtResponse).
        if (split && a.kind==K::HubbardLinearResponse && !handedOver && out.converged)
        {
            handedOver=true;
            bool real=false;
            for (size_t i=0;i<calc->Basis().GetNumIBS();i++) real = real || calc->Basis().GetRealIBS(i)!=nullptr;
            if (!real)
            {
                const std::string decided="ONE run serves the responses: the ground state has no real TRIM block to hand over";
                policy::Report("ResponseSplit", decided, policy::Verdict::Choose, "every k of this mesh is complex already");
                std::cout<<"[policy] "<<decided<<std::endl;
            }
            else
            {
                std::string path=spec.state.save;   // the run saved its final state there after its last Converge
                if (path.empty())
                {
                    fs::create_directories(outDir/"states");
                    path=(outDir/"states"/(stem+".groundstate.h5")).string();
                    calc->SaveState(path);
                }
                const std::string decided="TWO runs: the real ground state hands its state to a complex response run (no SCF in it)";
                const std::string reason="the ansatz is decided once per run from k and q: the ground state keeps its real TRIM blocks (faster, "
                                         "lighter), the response needs complex ones; run #1 is destroyed before run #2 is built, and run #2 "
                                         "adopts the state, verified by one Fock build ("+path+")";
                policy::Report("ResponseSplit", decided, policy::Verdict::Choose, reason);
                std::cout<<"[policy] "<<decided<<" -- "<<reason<<std::endl;
                calc.reset();                 // run #1's objects die with it
                qchem::report::Scope run2("responseRun");   // run #2's build and its adoption report one level in, under the hand-over
                calc=AdoptOf(spec, path);     // run #2
            }
        }
        switch (a.kind)
        {
        case K::EstimateHubbardU:
        {
            json est=json::array(); std::ostringstream os;
            for (const auto& e : calc->EstimateHubbardU())
            {
                est.push_back({{"site",e.site},{"l",e.l},{"Ubar_eV",e.Ubar*kHaToEV},{"Jbar_eV",e.Jbar*kHaToEV},{"Ueff_eV",e.Ueff()*kHaToEV}});
                os<<" site"<<e.site<<" l="<<e.l<<" Ueff="<<e.Ueff()*kHaToEV<<" eV;";
            }
            rj["estimates"]=est; pr.summary="ACBN0:"+os.str(); break;
        }
        case K::HubbardLoop:
        {
            SolidCalculation::HubbardLoop lp; lp.maxOuter=a.maxOuter; lp.tolU_eV=a.tolU_eV; lp.estimator=a.estimator;
            lp.nq=a.nq; lp.krylov=KrylovParams{.tol=a.tol}; lp.krylov.maxIter=a.maxIter; lp.krylov.restart=a.restart; lp.perturbed=a.perturb;
            auto L=calc->ConvergeHubbardU(stages.back().params, lp);
            rj["estimator"]=NameOf(a.estimator,kLoopEstimator);
            rj["outer"]=L.outer; rj["converged"]=L.converged; rj["scfConverged"]=L.scfConverged; rj["U_eV"]=L.U_eV;
            rj["estimated"]=L.estimated; rj["Uin_eV"]=L.Uin_eV; rj["moment"]=L.moment;
            json gap=json::array(); for (double g : L.gap_eV) gap.push_back(std::isfinite(g) ? json(g) : json(nullptr)); rj["gap_eV"]=gap;
            if (!L.failure.empty()) rj["failure"]=L.failure;
            pr.ok=L.converged && L.scfConverged && L.failure.empty();
            std::ostringstream os; os<<NameOf(a.estimator,kLoopEstimator)<<": "<<L.outer<<" outer steps, U "<<(L.converged?"CONVERGED":"NOT converged")
                                  <<", last SCF "<<(L.scfConverged?"converged":"NOT converged");
            if (!L.failure.empty()) os<<"; FAILED: "<<L.failure;
            pr.summary=os.str();
            if (auto R=calc->Result()) { out.converged=true; out.energy=calc->LastIterateTerms().GetTotalEnergy(); }   // the final-U SCF is now the run's answer
            break;
        }
        case K::IndependentResponse:
        {
            auto R=calc->IndependentResponse(a.nq);
            pr.ok=bool(R);
            if (R) { rj["labels"]=R->labels; rj["gap"]=R->gap; pr.summary="chi0 over "+std::to_string(R->labels.size())+" channels, gap "+std::to_string(R->gap); }
            else   { pr.summary="FAILED: "+R.Error().detail; rj["failure"]=R.Error().detail; }
            break;
        }
        case K::HubbardLinearResponse:
        {
            KrylovParams kp{.tol=a.tol}; kp.maxIter=a.maxIter; kp.restart=a.restart;
            auto R=calc->HubbardLinearResponse(a.nq, kp, a.perturb, a.resolve, a.qStars, a.littleGroupK);
            pr.ok=bool(R);
            if (R)
            {
                // Per q: chi0 and chi (measured x perturbed; real and imaginary parts -- off-diagonal chi(q) is complex at a non-TRIM
                // q), the Krylov residuals and costs.  Then the supercell's real-space perturbed blocks and U, hp.x's quantities.
                const auto& s0=R->atQ[0];
                rj["nq"]=QMeshJson(a.nq); rj["labels"]=s0.labels; rj["perturbed"]=s0.perturbed;
                auto mat=[](const auto& M, auto part)
                { json m=json::array(); for (size_t i=0;i<M.rows();++i) { json row=json::array(); for (size_t c=0;c<M.columns();++c) row.push_back(part(M(i,c))); m.push_back(row); } return m; };
                auto re=[](const auto& z){ return z.real(); };
                auto im=[](const auto& z){ return z.imag(); };
                json perQ=json::array();
                for (size_t iq=0;iq<R->atQ.size();++iq)
                {
                    const auto& s=R->atQ[iq]; const rvec3_t q=R->q[iq].q();
                    perQ.push_back({{"q",{q.x,q.y,q.z}},{"gap",s.gap},{"residual",s.residual},{"iterations",s.iterations},
                                    {"chi0_real",mat(s.chi0,re)},{"chi0_imag",mat(s.chi0,im)},{"chi_real",mat(s.chi,re)},{"chi_imag",mat(s.chi,im)}});
                }
                rj["perQ"]=perQ;
                rj["chi0R"]=mat(R->Chi0RealSpace(),[](double x){return x;});
                rj["chiR"] =mat(R->ChiRealSpace (),[](double x){return x;});
                const rvec_t U=R->U(); std::vector<double> UeV; for (size_t c=0;c<U.size();++c) UeV.push_back(ToEV(U[c]));
                rj["U_eV"]=UeV;
                if (U.size()>1)                                   // the orbital-resolved 1x1 inversion (Macke 2024's pointwise) beside hp.x's
                {
                    const rvec_t us=R->UScalar(); std::vector<double> u1; for (size_t c=0;c<us.size();++c) u1.push_back(ToEV(us[c]));
                    rj["U1x1_eV"]=u1;
                    if (!R->group.empty()) rj["manifold"]=R->group;
                }
                json V=json::array();
                for (const auto& v : calc->InterSiteHubbardV(*R))
                    V.push_back({{"from",v.from},{"to",v.to},{"fromLabel",s0.labels[v.from]},{"toLabel",s0.labels[v.to]},
                                 {"R",{v.R.x,v.R.y,v.R.z}},{"d_bohr",v.distance},{"V_eV",ToEV(v.V)}});
                rj["V"]=V;
                std::ostringstream os; os<<"self-consistent chi over "<<s0.labels.size()<<" channels at "<<R->atQ.size()<<" q; U(eV)";
                for (double u : UeV) os<<" "<<u;
                pr.summary=os.str();
            }
            else { pr.summary="FAILED: "+R.Error().detail; rj["failure"]=R.Error().detail; }
            break;
        }
        case K::HubbardFiniteDifference:
        {
            std::vector<size_t> J=a.perturb; if (J.empty()) J.push_back(0);
            json fds=json::array(); bool ok=true;
            for (size_t j : J)
            {
                auto R=calc->HubbardFiniteDifferenceChi(j, a.alpha, stages.back().params);
                if (R) { std::vector<double> c; for (size_t i=0;i<R->chi.size();++i) c.push_back(R->chi[i]); fds.push_back({{"perturbed",j},{"chi",c},{"restored",R->restored}}); ok = ok && R->restored; }
                else   { fds.push_back({{"perturbed",j},{"failure",R.Error().details}}); ok=false; }
            }
            rj["fd"]=fds; pr.ok=ok; pr.summary=std::to_string(J.size())+" manifold(s) perturbed"; break;
        }
        case K::Displace:
        {
            json pj=json::array(); std::ostringstream os;
            const auto D=Displace(spec, a.points);
            pr.ok=D.ok;
            rj["reference"]={{"energy",D.E0 ? json(*D.E0) : json(nullptr)},{"iterations",D.it0}};
            if (!D.failure0.empty()) rj["reference"]["failure"]=D.failure0;
            for (size_t k=0;k<D.points.size();++k)
            {
                const auto& q=D.points[k];
                json x={{"atoms",a.points[k].atoms},{"delta",{a.points[k].delta.x,a.points[k].delta.y,a.points[k].delta.z}},{"iterations",q.iterations}};
                if (q.energy) { x["energy"]=*q.energy; if (D.E0) x["dE"]=*q.energy-*D.E0; }
                else          x["failure"]=q.failure;
                pj.push_back(x);
            }
            rj["points"]=pj;
            if (D.E0) os<<"reference "<<std::setprecision(12)<<*D.E0<<" Ha;"; else os<<"reference FAILED: "<<D.failure0<<";";
            os<<std::setprecision(6);
            for (const auto& q : D.points) { if (q.energy && D.E0) os<<" dE "<<std::scientific<<(*q.energy-*D.E0); else os<<" FAILED"; }
            pr.summary=os.str(); break;
        }
        case K::GridAudit:
        case K::EggBox:
        case K::XCMeshAudit:
        {
            auto xcRef=[&]
            {
                const auto S=calc->SizedSetup();
                return qcMesh::XCAuditReference(S->xcMesh, a.xcRefNRadial, a.xcRefAngularDegree, S->opts.solid.beckeLadder);
            };
            const auto O= a.kind==K::EggBox      ? calc->EggBox(a.shift)
                        : a.kind==K::XCMeshAudit ? calc->XCMeshAudit(xcRef())
                        :                          calc->GridAudit(a.referenceMargin);
            if (!O)
            {   // REFUSED ON RESOURCES (the audit RSS policy): nothing was built; the record carries the numbers and the escape
                const ResourceRefusal& r=O.Error();
                rj["refused"]=Describe(r.footprint, r.available); rj["escape"]=r.escape;
                pr.ok=false; pr.summary="REFUSED ON RESOURCES: "+r.Summary(); break;
            }
            const auto& R=O.Value();
            json t=json::array();
            for (const auto& x : R.terms) t.push_back({{"term",x.label},{"entries",x.entries},{"maxAbsDH",x.maxAbsΔH},{"dE",x.ΔE},
                                                       {"x",x.x},{"bound",x.bound},{"eta",x.eta},{"etaOn",x.etaOn},
                                                       {"trDdH",x.trDdH},{"sNormDH",x.sNormDH}});
            rj["reference"]=R.reference; rj["terms"]=t; rj["maxAbsDH"]=R.maxAbsΔH; rj["dE"]=R.ΔE;
            auto vj=[](const auto& vs) { json o=json::object(); for (const auto& v : vs) o[v.audit]={{"status",v.status},{"alias",v.alias}}; return o; };
            rj["verdicts"]=vj(R.runVerdicts); rj["referenceVerdicts"]=vj(R.referenceVerdicts);
            // ok = the audit MEETS THE RUN'S OWN TARGET on the run's own measure (lane H beta (c), 2026-11-07): max|dH_ij| <= eps for a client
            // that consumes matrix elements, |dE| <= eps N for an energy run -- the pilot's Target().  It used to stay the default true,
            // whatever the audit measured.  A run with no target (accuracy 0, stated rasters) has nothing to judge against: ok, and says so.
            const auto S=calc->SizedSetup();
            const ClientConsumes c=calc->Consumes();
            const double target=S->opts.accuracy.Target(c, S->structure->GetNumAtoms());
            const bool onH = c==ClientConsumes::MatrixElements;
            const double measured = onH ? R.maxAbsΔH : std::abs(R.ΔE);
            pr.ok = target<=0.0 || measured<=target;
            rj["target"]=target>0.0 ? json(target) : json(nullptr); rj["measure"]=onH ? "maxAbsDH" : "absDE";
            std::ostringstream os; os<<"vs "<<R.reference<<": max|dH_ij| "<<R.maxAbsΔH<<", dE "<<R.ΔE<<" Ha";
            if (target>0.0) os<<(pr.ok ? "  -- meets " : "  -- OVER ")<<(onH ? "max|dH_ij| " : "|dE| ")<<target<<" Ha (the run's target)";
            else            os<<"  (no accuracy target: not judged)";
            pr.summary=os.str(); break;
        }
        }
        rj["ok"]=pr.ok; rj["summary"]=pr.summary; results.push_back(rj);
        out.postSCF.push_back(pr);
    }
    deckTiming.Report();   // the run's `timing`, at its root: sorted by cost + PEAK RSS + SETUP / SCF TOTAL (what scripts/retake5a greps)
    if (!results.empty() || refused)
    {   // fold the results into the record (the revision written before the run keeps the deck; `results` is what the run produced)
        std::ifstream in(out.revision); json rec=json::parse(in); in.close();
        rec["results"]={{"converged",out.converged},{"postSCF",results}};
        if (refused) rec["results"]["refusedOnResources"]=Describe(refused->footprint, refused->available);
        if (out.energy) rec["results"]["energy"]=*out.energy;
        std::ofstream os(out.revision,std::ios::trunc); os<<rec.dump(2)<<"\n";
    }
    return out;
}

//=== Build: the calculation a deck states, stopped before its SCF =====================================
std::vector<SCFStage> ScheduleToRun(const RunSpec& spec)
{
    std::vector<SCFStage> stages;
    for (const auto& s : spec.Stages()) stages.push_back({s.scf,s.accelerator,s.recipe});
    if (!spec.state.restartFrom.empty()) stages={stages.back()};   // a restart is ONE stage: the schedule's final one
    return stages;
}

// THE SETUP, resolved ONCE: the basis from the deck's recipe (the trim stated, or vetted at the run's own orthoTol and frozen), and every
// decision the run makes before it builds anything.  The first stage's SCF settings, accelerator and recipe (a restart's: its one stage's)
// are the setup's, so Prepare builds exactly stage 0; each later stage names its own.  \a spec is Build's (its family copied into solid, its save path resolved).
// A deck is USER INPUT, so its setup goes through SolidSetup::Admit (OpenWork "SolidSetup refusals are throws", 2026-11-14): Run reports a
// refusal as its outcome and files it in the record; the doors whose signatures carry no refusal (Build, View, MakeBasisShells: the GUI
// binding calls them) say it as the deck's own refusal, a throw naming the policy -- the deck's convention for every key it refuses.
static Outcome<std::shared_ptr<const SolidSetup>, SetupRefusal> AdmitSetup(const RunSpec& spec, const Materials::Material& mat,
                                                                           const std::vector<SCFStage>& stages, RunPart part=RunPart::Whole)
{
    Lattice_3D lat(*mat.cell, spec.kmesh);
    RunOptions so=spec.opts; so.scf=stages.front().params; so.accelerator=stages.front().accelerator; so.recipe=stages.front().recipe;
    so.solid.responseQMeshes = part==RunPart::GroundState ? std::vector<ivec3_t>{} : PlannedResponseQMeshes(spec);
    return SolidSetup::Admit(lat, RecipeOf(spec.basis), so, spec.accel);
}
static std::string Refusal(const SetupRefusal& r) {return "deck: the setup refused the deck (policy "+r.policy+"): "+r.reason;}
static std::string Refusal(const BuildRefusal& r) {return "deck: the build refused the deck (policy "+r.policy+"): "+r.reason;}
static std::shared_ptr<const SolidSetup> SetupOf(const RunSpec& spec, const Materials::Material& mat, const std::vector<SCFStage>& stages,
                                                 RunPart part)
{
    auto S=AdmitSetup(spec, mat, stages, part);
    if (!S) throw std::runtime_error(Refusal(S.Error()));
    return S.TakeValue();
}

std::unique_ptr<SolidCalculation> Build(const RunSpec& specIn, SolidCalculation::Observer onIteration)
{
    auto C=PrepareOf(specIn, std::move(onIteration));
    if (!C) throw std::runtime_error(Refusal(C.Error()));
    return C.TakeValue();
}

std::unique_ptr<Calculation> Prepare(const RunSpec& spec)
{
    switch (KindOf(spec))
    {
    case RunKind::Solid: return Build(spec);
    case RunKind::Molecule:
    {
        auto setup=MoleculeSetup::Admit(StructureData::GetMolecule(spec.structure), spec.opts, spec.accel);
        if (!setup) throw std::runtime_error(Refusal(setup.Error()));
        auto built=MoleculeCalculation::Prepare(setup.Value());
        if (!built) throw std::runtime_error(Refusal(built.Error()));
        return built.TakeValue();
    }
    case RunKind::Atom:
    {
        auto setup=AtomSetup::Admit(StructureData::GetMolecule(spec.structure), spec.opts, spec.accel);   // ONE atom: its Z and charge are the structure's
        if (!setup) throw std::runtime_error(Refusal(setup.Error()));
        auto built=AtomCalculation::Prepare(setup.Value());
        if (!built) throw std::runtime_error(Refusal(built.Error()));
        return built.TakeValue();
    }
    }
    throw std::logic_error("deck::Prepare: RunKind");
}

static Outcome<std::unique_ptr<SolidCalculation>, SetupRefusal> PrepareOf(const RunSpec& specIn, SolidCalculation::Observer onIteration,
                                                                          RunPart part)
{
    using R=Outcome<std::unique_ptr<SolidCalculation>, SetupRefusal>;
    namespace fs=std::filesystem;
    RunSpec spec=specIn;
    const Materials::Material mat=Resolve(spec);
    const std::string& r=spec.state.restartFrom;
    if (!r.empty())
    {
        const bool isPath = r.find('/')!=std::string::npos || (r.size()>3 && r.substr(r.size()-3)==".h5");
        if (!isPath) throw std::runtime_error("deck: state.restartFrom '"+r+"' is a revision STEM, which only Run (with its outDir) can resolve; "
                                              "Build takes the saved state's path");
        if (!fs::exists(r)) throw std::runtime_error("deck: state.restartFrom '"+r+"': no saved state there");
    }
    spec.opts.solid.orbitalBasis=spec.basis.family; spec.opts.solid.planeWaveEcut=spec.basis.ecut;   // the basis block names the family; the facade reads it from solid
    if (!spec.state.save.empty())
    {
        if (spec.state.save=="auto") throw std::runtime_error("deck: state.save 'auto' names a file under Run's outDir; Build takes a path");
        spec.opts.saveStateTo=spec.state.save;
    }
    const std::vector<SCFStage> stages=ScheduleToRun(spec);
    auto admitted=AdmitSetup(spec, mat, stages, part);
    if (!admitted) return R::Fail(admitted.Error());
    const auto setup=admitted.TakeValue();
    if (r.empty())
    {   // a BUILD refusal is reported as the setup's is (Run files both alike: a policy row and its sentence, no SCF ran)
        auto built=SolidCalculation::Prepare(setup, std::move(onIteration));
        if (!built) return R::Fail(SetupRefusal{built.Error().policy, built.Error().reason});
        return R::Ok(built.TakeValue());
    }
    // a refusal THROWS -- a deck never silently re-seeds
    auto P=SolidCalculation::PrepareRestart(r, setup, stages.front().params, std::move(onIteration));
    if (!P) throw std::runtime_error("deck: state.restartFrom "+r+" refused: "+P.Error().details);
    return R::Ok(P.TakeValue());
}

static std::unique_ptr<SolidCalculation> AdoptOf(const RunSpec& specIn, const std::string& statePath)
{
    RunSpec spec=specIn;
    const Materials::Material mat=Resolve(spec);
    spec.opts.solid.orbitalBasis=spec.basis.family; spec.opts.solid.planeWaveEcut=spec.basis.ecut;   // as Build
    spec.opts.saveStateTo.clear();   // run #2 converges nothing to save (a fallback restart re-saves nothing either)
    const std::vector<SCFStage> stages=ScheduleToRun(spec);
    auto A=SolidCalculation::Adopt(statePath, SetupOf(spec, mat, stages, RunPart::Response), stages.back().params);
    if (!A) throw std::runtime_error("deck: the response run could not adopt the ground state "+statePath+": "+A.Error().details);
    return A.TakeValue();
}

Outcome<std::unique_ptr<StateView>, RestartRefusal> View(const RunSpec& specIn, const std::string& statePath)
{
    RunSpec spec=specIn;
    const Materials::Material mat=Resolve(spec);
    spec.opts.solid.orbitalBasis=spec.basis.family; spec.opts.solid.planeWaveEcut=spec.basis.ecut;   // as Build
    spec.opts.saveStateTo.clear();   // a view saves nothing
    return SolidCalculation::ViewState(statePath, SetupOf(spec, mat, ScheduleToRun(spec)));
}

//=== The shell what-if: the deck's setup, its UNTRIMMED basis, and the shells its trim keeps =======================================
BasisShells MakeBasisShells(RunSpec spec)
{
    const Materials::Material mat=Resolve(spec);
    spec.opts.solid.orbitalBasis=spec.basis.family; spec.opts.solid.planeWaveEcut=spec.basis.ecut;   // as Build: the basis block names the family
    BasisShells out;
    out.setup=SetupOf(spec, mat, ScheduleToRun(spec));
    out.table=out.setup->MakeShellConditioning();
    out.floor=spec.opts.orthoTol>0.0 ? spec.opts.orthoTol : ConditioningFloor;
    out.deckKept=out.table->KeptBy(out.setup->recipe->trim);              // the trim FROZEN in the setup: stated, or autoTrim's choice
    out.autoTrimKept=out.table->KeptBy(out.setup->MakeVetTrim(out.floor));
    return out;
}

BasisSet::Lattice::ShellConditioning::Answer BasisShells::Analyze(std::vector<bool> kept) const
{
    if (kept.empty()) kept=deckKept;
    auto a=table->Analyze(kept, floor);
    report::Set("deckKept", deckKept);
    report::Set("autoTrimKept", autoTrimKept);
    return a;
}

} // namespace qchem::deck
