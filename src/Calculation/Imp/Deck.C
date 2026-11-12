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
import qchem.PeriodicTable;                 // the material's element symbols, for an inlined basis
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94Import;   // SubsetBasis: the elements a deck inlines
import qchem.Mesh.XCPolicy;                 // ResolveXCMesh (a recipe field under cellKind Auto is refused at resolve)
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.StructureData;   // KindOf

namespace qchem::deck
{
namespace
{

//! Reads ONE JSON object into a struct: each \c Get is optional (a missing key keeps the default), and \c Done REJECTS any key
//! no \c Get asked for -- the typo guard.  \a path names the object in the message (`solid.xcMesh`).
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
const auto kAngularForm = NAMES(BasisSet::Gaussian::Angular, {"cartesian",BasisSet::Gaussian::Angular::Cartesian}, {"spherical",BasisSet::Gaussian::Angular::Spherical});
const auto kFamily = NAMES(SolidCalcOptions::OrbitalBasis, {"gaussian",SolidCalcOptions::OrbitalBasis::GPW}, {"planeWave",SolidCalcOptions::OrbitalBasis::PlaneWave});
const auto kVxc    = NAMES(Hamiltonian::VxcFit, {"Auto",Hamiltonian::VxcFit::Auto}, {"PlaneWave",Hamiltonian::VxcFit::PlaneWave}, {"Delta",Hamiltonian::VxcFit::Delta}, {"Auxiliary",Hamiltonian::VxcFit::Auxiliary});
const auto kRho    = NAMES(Hamiltonian::XCρ, {"exact",Hamiltonian::XCρ::Exact}, {"fit",Hamiltonian::XCρ::Fit});
const auto kAcc    = NAMES(AccType, {"DIIS",AccType::DIIS}, {"GDM",AccType::GDM}, {"Ladder",AccType::Ladder}, {"Null",AccType::Null},
                         {"Auto",AccType::Auto});
const auto kSignal = NAMES(SCFAccelerators::ScheduleSignal, {"Error",SCFAccelerators::ScheduleSignal::Error}, {"EnergyChange",SCFAccelerators::ScheduleSignal::EnergyChange});
const auto kSeed   = NAMES(SeedStrategy, {"Default",SeedStrategy::Default}, {"CoreGuess",SeedStrategy::CoreGuess}, {"Uniform",SeedStrategy::Uniform},
                                         {"SAD",SeedStrategy::SAD}, {"IonicSAD",SeedStrategy::IonicSAD});
const auto kOrtho  = NAMES(Ortho, {"Cholesky",Cholesky}, {"Eigen",Eigen}, {"SVD",SVD}, {"Auto",Auto}, {"CholeskyPivoted",CholeskyPivoted}, {"Orthonormal",Orthonormal});
const auto kCell   = NAMES(UnitCellKind, {"Uniform",UnitCellKind::Uniform}, {"Becke",UnitCellKind::Becke}, {"Auto",UnitCellKind::Auto});
const auto kResolution = NAMES(Response::HubbardResolution, {"manifold",Response::HubbardResolution::Manifold},
                               {"slot",Response::HubbardResolution::Slot}, {"eigen",Response::HubbardResolution::Eigen});
const auto kLoopEstimator = NAMES(SolidCalculation::HubbardLoop::Estimator, {"ACBN0",SolidCalculation::HubbardLoop::Estimator::ACBN0},
                                  {"linearResponse",SolidCalculation::HubbardLoop::Estimator::LinearResponse});
const auto kProjector = NAMES(Hamiltonian::HubbardProjectorKind, {"siteShells",Hamiltonian::HubbardProjectorKind::SiteShells},
                              {"pseudoAtomic",Hamiltonian::HubbardProjectorKind::PseudoAtomic});
const auto kRadial = NAMES(RadialKind, {"MHL",RadialKind::MHL}, {"Log",RadialKind::Log}, {"Linear",RadialKind::Linear});
const auto kAngular= NAMES(AngularKind, {"Lebedev",AngularKind::Lebedev}, {"GaussLegendre",AngularKind::GaussLegendre});
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
            {"rasterX",t.rasterX},{"vlocRasterX",t.vlocRasterX},{"xcRasterX",t.xcRasterX}};
}
void FromJson(const json& j, BasisSet::Gaussian::GPWTolerances& t)
{
    Fields f(j,"tolerances");
    f.Get("vlocEps",t.vlocEps); f.Get("localPPRelCutoff",t.localPPRelCutoff); f.Get("relFieldSharp",t.relFieldSharp);
    f.Get("mgridEcuts",t.mgridEcuts); f.Get("screenEps",t.screenEps); f.Get("fieldSharp",t.fieldSharp);
    f.Get("relCutoff",t.relCutoff); f.Get("densityEps",t.densityEps);
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
    f.Done();
}

//=== SolidAcceleratorOptions ===================================================================================================
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
json ToJson(const SCFAccelerators::SolidAcceleratorOptions& o)
{
    return {{"diis",DiisJson(o.diis)},{"gdm",GdmJson(o.gdm)},
            {"ladder",{{"gdm",GdmJson(o.ladderGdm)},{"eThresh",o.ladderEThresh},{"stall",o.ladderStall},{"floor",o.ladderFloor},
                       {"switchAt",o.ladderSwitchAt},{"signal",NameOf(o.ladderSignal,kSignal)}}}};
}
void FromJson(const json& j, SCFAccelerators::SolidAcceleratorOptions& o)
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

//=== MeshParams =================================================================================
json ToJson(const qcMesh::MeshParams& m)
{
    return {{"radial",NameOf(m.radial,kRadial)},{"nRadial",m.nRadial},{"mhlM",m.mhl_m},{"mhlAlpha",m.mhl_alpha},{"logStart",m.logStart},
            {"logStop",m.logStop},{"angular",NameOf(m.angular,kAngular)},{"angularDegree",m.angularDegree},{"angRot",m.angRot},
            {"beckeOrder",m.beckeOrder},{"nUniform",m.nUniform},{"eCut",m.eCut},{"relCutoff",m.relCutoff},
            {"cellKind",NameOf(qcMesh::UnitCellKind(m.cellKind),kCell)},{"beckeEps",m.beckeEps}};
}
void FromJson(const json& j, qcMesh::MeshParams& m)
{
    Fields f(j,"xcMesh");
    f.GetEnum("radial",m.radial,kRadial); f.Get("nRadial",m.nRadial); f.Get("mhlM",m.mhl_m); f.Get("mhlAlpha",m.mhl_alpha);
    f.Get("logStart",m.logStart); f.Get("logStop",m.logStop); f.GetEnum("angular",m.angular,kAngular);
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

//=== SolidCalcOptions ===========================================================================
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

json ToJson(const SolidCalcOptions& o)
{
    json species=json::array();
    for (const auto& s : o.species) species.push_back({s.first,s.second});
    json hub=json::array();
    for (const auto& h : o.hubbard) hub.push_back(ToJson(h));
    json out={{"Nelec",o.Nelec},{"charge",o.charge},{"multiplicity",o.multiplicity},{"pp",o.pp},{"species",species},
            {"densityEcut",o.densityEcut},{"cutoffFactor",o.cutoffFactor},{"ladderFactor",o.ladderFactor},
            {"images",NameOf(o.images,kImages)},{"kShift",{o.kShift.x,o.kShift.y,o.kShift.z}},
            {"xcMesh",ToJson(o.xcMesh)},{"tolerances",ToJson(o.tolerances)},
            {"accuracy",o.accuracy.eps},{"policy",ToJson(o.policy)},{"hubbard",hub},
            {"accelerator",NameOf(o.accelerator,kAcc)},{"perKFermiLevels",o.perKFermiLevels},{"imposeSymmetry",o.imposeSymmetry},{"xcMeshSymmetric",o.xcMeshSymmetric},
            {"seed",NameOf(o.seed,kSeed)},{"ortho",NameOf(o.ortho,kOrtho)},{"orthoTol",o.orthoTol},{"numBands",o.numBands},{"forceComplex",o.forceComplex},
            {"spinsShareFermi",o.spinsShareFermi},{"greyImposition",o.greyImposition},{"momFromSeed",o.momFromSeed},
            {"siteSpins",o.siteSpins},{"siteCharges",o.siteCharges},{"label",o.label},{"saveStateTo",o.saveStateTo}};
    if (o.accuracy.eta) out["accuracyCalibration"]={{"vlocLong",o.accuracy.eta->vlocLong},{"density",o.accuracy.eta->density}};   // stated only
    if (o.accuracy.vlocRoutingEta!=BasisSet::Gaussian::GridAccuracy{}.vlocRoutingEta) out["accuracyCalibration"]["vlocRouting"]=o.accuracy.vlocRoutingEta;
    if (o.accuracy.rungGateC!=BasisSet::Gaussian::GridAccuracy{}.rungGateC) out["accuracyCalibration"]["rungGateC"]=o.accuracy.rungGateC;
    if (o.accuracy.relCutoffSafety!=BasisSet::Gaussian::GridAccuracy{}.relCutoffSafety) out["accuracyCalibration"]["relCutoffSafety"]=o.accuracy.relCutoffSafety;
    if (o.accuracy.samplerEta!=BasisSet::Gaussian::GridAccuracy{}.samplerEta) out["accuracyCalibration"]["samplerEta"]=o.accuracy.samplerEta;
    if (o.accuracy.samplerResidualMargin!=BasisSet::Gaussian::GridAccuracy{}.samplerResidualMargin)
        out["accuracyCalibration"]["samplerResidualMargin"]=o.accuracy.samplerResidualMargin;
    if (o.auditFootprint!=kAuditFootprintCalibration)
        out["auditFootprint"]={{"eta",o.auditFootprint.eta},{"kappa",o.auditFootprint.kappa},{"bandLow",o.auditFootprint.bandLow},
                               {"bandHigh",o.auditFootprint.bandHigh}};
    if (o.beckeLadder!=qcMesh::BeckeLadderCalibration{})
        out["beckeLadder"]={{"angularError53",o.beckeLadder.angularError53},{"angularFactor",o.beckeLadder.angularFactor},
                            {"radialError100",o.beckeLadder.radialError100},{"radialFactor",o.beckeLadder.radialFactor},
                            {"radialShare",o.beckeLadder.radialShare},{"angRot",o.beckeLadder.angRot}};
    return out;
}
void FromJson(const json& j, SolidCalcOptions& o)
{
    Fields f(j,"solid");
    f.Get("Nelec",o.Nelec); f.Get("charge",o.charge); f.Get("multiplicity",o.multiplicity); f.Get("pp",o.pp);
    if (j.contains("species"))
    {
        o.species.clear();
        for (const auto& s : f.At("species"))
        {
            if (!s.is_array() || s.size()!=2) throw std::runtime_error("deck: 'solid.species' entries are [\"Element\", valence], got "+s.dump());
            o.species.emplace_back(s.at(0).get<std::string>(), s.at(1).get<int>());
        }
    }
    else f.At("species");
    f.Get("densityEcut",o.densityEcut); f.Get("cutoffFactor",o.cutoffFactor); f.Get("ladderFactor",o.ladderFactor);
    if (j.contains("raster"))
        throw std::runtime_error("deck: 'solid.raster' is retired (C2 step 3, doc/GridEpsilonPlan.md §1): the rasters are SIZED from "
                                 "'solid.accuracy' -- drop the key; to give them yourself state \"accuracy\": 0 and tolerances.rasterX "
                                 "(0 = the ball's own raster, the old BallOnly)");
    f.GetEnum("images",o.images,kImages);
    if (j.contains("kShift"))
    {
        const json& k=f.At("kShift");
        if (!k.is_array() || k.size()!=3) throw std::runtime_error("deck: 'solid.kShift' must be [x,y,z], got "+k.dump());
        o.kShift=rvec3_t(k.at(0).get<double>(),k.at(1).get<double>(),k.at(2).get<double>());
    }
    else f.At("kShift");
    if (j.contains("xcMesh"))     FromJson(f.At("xcMesh"),o.xcMesh);         else f.At("xcMesh");
    if (j.contains("tolerances")) FromJson(f.At("tolerances"),o.tolerances); else f.At("tolerances");
    f.Get("accuracy",o.accuracy.eps);
    if (!(o.accuracy.eps>=0.0)) throw std::runtime_error("deck: 'solid.accuracy' (Ha per atom on an energy run, per matrix element on a geometry run) "
                                                         "must be >= 0 (0 = no target: the rasters tolerances.rasterX ... state)");
    if (j.contains("beckeLadder"))   // STATED (expert): where a geometry run's pilot enters the Becke ladder
    {
        Fields c(f.At("beckeLadder"),"solid.beckeLadder");
        qcMesh::BeckeLadderCalibration& b=o.beckeLadder;
        c.Get("angularError53",b.angularError53); c.Get("angularFactor",b.angularFactor); c.Get("radialError100",b.radialError100);
        c.Get("radialFactor",b.radialFactor); c.Get("radialShare",b.radialShare); c.Get("angRot",b.angRot);
        c.Done();
        if (!(b.angularError53>0.0 && b.radialError100>0.0))
            throw std::runtime_error("deck: 'solid.beckeLadder' errors are measured max|dH_xc| values, > 0");
        if (!(b.angularFactor>0.0 && b.angularFactor<1.0 && b.radialFactor>0.0 && b.radialFactor<1.0 && b.radialShare>0.0 && b.radialShare<=1.0))
            throw std::runtime_error("deck: 'solid.beckeLadder' factors are falls per rung, in (0, 1), and radialShare is in (0, 1]");
    }
    if (j.contains("auditFootprint"))   // STATED (expert): the audit RSS policy's calibration (absent = its measured record)
    {
        Fields c(f.At("auditFootprint"),"solid.auditFootprint");
        AuditFootprintCalibration& a=o.auditFootprint;
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
        c.Get("vlocRouting",o.accuracy.vlocRoutingEta);
        c.Get("rungGateC",o.accuracy.rungGateC);
        c.Get("relCutoffSafety",o.accuracy.relCutoffSafety);
        c.Get("samplerEta",o.accuracy.samplerEta);
        c.Get("samplerResidualMargin",o.accuracy.samplerResidualMargin);
        c.Done();
        if (!(o.accuracy.rungGateC>=0.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.rungGateC' is the C of the ladder's rung gate, >= 0 (0 = no completion rung)");
        if (!(o.accuracy.relCutoffSafety>0.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.relCutoffSafety' is the stiffness of the pair->level rule, > 0 (default 2)");
        if (!(o.accuracy.samplerEta>0.0 && o.accuracy.samplerEta<=1.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.samplerEta' is the rho-tilde NUFFT's delivered/bound, in (0, 1] (1 = rigorous)");
        if (!(o.accuracy.samplerResidualMargin>0.0 && o.accuracy.samplerResidualMargin<=1.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.samplerResidualMargin' is the share of the SCF density gate the sampling may use, in (0, 1]");
        if (!(eta.vlocLong>0.0 && eta.vlocLong<=1.0 && eta.density>0.0 && eta.density<=1.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration' factors are eta = delivered/bound, in (0, 1] (1 = rigorous)");
        if (!(o.accuracy.vlocRoutingEta>0.0))
            throw std::runtime_error("deck: 'solid.accuracyCalibration.vlocRouting' is V_loc long's routing calibration eta = delivered/bound, > 0 "
                                     "(a stated 'solid.tolerances.vlocEps' wins over it)");
        if (etaStated) o.accuracy.eta=eta;
    }
    if (j.contains("policy"))     FromJson(f.At("policy"),o.policy);         else f.At("policy");
    if (j.contains("vxcFit"))             // RETIRED 2026-11-10 (V-XC-ONE-TERM step 3): the fit basis is one of the run's XC choices
        throw std::runtime_error("deck: 'solid.vxcFit' is retired (2026-11-10): the Vxc fit basis moved into the run's 'xc' object -- "
                                 "state \"xc\": {\"fit\": \""+(j.at("vxcFit").is_string() ? j.at("vxcFit").get<std::string>() : std::string("Delta"))
                                 +"\"} (beside \"rho\": \"exact\" | \"fit\")");
    if (j.contains("hubbard"))
    {
        o.hubbard.clear();
        size_t i=0;
        for (const auto& h : f.At("hubbard")) { o.hubbard.emplace_back(); FromJson(h,o.hubbard.back(),"solid.hubbard."+std::to_string(i++)); }
    }
    else f.At("hubbard");
    if (j.contains("globalFermi"))           // RETIRED 2026-10-21 (D-FERMI-DEFAULT): never silently re-read with the opposite default
        throw std::runtime_error("deck: solid.globalFermi is RETIRED -- one Fermi level across the k-mesh is now the DEFAULT for a "
                                 "smeared run; delete the key (true), or set solid.perKFermiLevels true for the old per-k default (false)");
    f.GetEnum("accelerator",o.accelerator,kAcc); f.Get("perKFermiLevels",o.perKFermiLevels); f.Get("imposeSymmetry",o.imposeSymmetry);
    f.Get("xcMeshSymmetric",o.xcMeshSymmetric);   // false = the plain atom-centred mesh: a displacement path needs it (SolidCalcOptions)
    f.GetEnum("seed",o.seed,kSeed); f.GetEnum("ortho",o.ortho,kOrtho); f.Get("orthoTol",o.orthoTol); f.Get("numBands",o.numBands); f.Get("forceComplex",o.forceComplex);
    f.Get("spinsShareFermi",o.spinsShareFermi); f.Get("greyImposition",o.greyImposition); f.Get("momFromSeed",o.momFromSeed);
    f.Get("siteSpins",o.siteSpins); f.Get("siteCharges",o.siteCharges); f.Get("label",o.label); f.Get("saveStateTo",o.saveStateTo);
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
    if (b.family!=SolidCalcOptions::OrbitalBasis::GPW) { j["family"]=NameOf(b.family,kFamily); j["ecut"]=b.ecut; }
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
                                 "\"autoTrim\": true beside a solid.orthoTol > 0, and otherwise drop the key");
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

json ToJson(const RunSpec& r)
{
    json j={{"structure", r.inlineStructure.is_null() ? json(r.structure) : r.inlineStructure}};
    j["kmesh"]={r.kmesh.x,r.kmesh.y,r.kmesh.z}; j["basis"]=BasisJson(r.basis); j["solid"]=ToJson(r.solid);
    // The ONE accuracy target is the RUN's (V-FACADE-ACCURACY): the record states it at the top level, never as `solid.accuracy`.
    j["accuracy"]=j["solid"]["accuracy"]; j["solid"].erase("accuracy");
    // THE XC CHOICES are the run's (like `accuracy`): which fit basis represents v_xc, and which density XC sees (user ruling 2026-11-09).
    j["xc"]={{"fit",NameOf(r.solid.vxcFit,kVxc)},{"rho",NameOf(r.solid.xcRho,kRho)}};
    j["acceleratorOptions"]=ToJson(r.accel);
    if (!r.state.save.empty() || !r.state.restartFrom.empty())
    {
        json st=json::object();
        if (!r.state.save.empty())        st["save"]=r.state.save;
        if (!r.state.restartFrom.empty()) st["restartFrom"]=r.state.restartFrom;
        j["state"]=st;
    }
    if (!r.postSCF.empty()) { json a=json::array(); for (const auto& x : r.postSCF) a.push_back(PostJson(x)); j["postSCF"]=a; }
    if (r.threads!=ThreadCounts{}) j["threads"]={{"openmp",r.threads.openmp},{"blas",r.threads.blas}};   // stated only: absent = the process's
    if (r.schedule.empty()) j["scf"]=ToJson(r.scf);                      // one stage: the record says it once
    else
    {
        json st=json::array();
        for (const auto& s : r.schedule) st.push_back({{"accelerator",NameOf(s.accelerator,kAcc)},{"scf",ToJson(s.scf)}});
        j["schedule"]=st;
    }
    return j;
}
void FromJson(const json& j, RunSpec& r)
{
    Fields f(j,"run");
    if (!j.contains("structure") || !(j.at("structure").is_string() || j.at("structure").is_object()))
        throw std::runtime_error("deck: 'structure' is required: a NAME from materials.json or molecules.json, or the cell itself "
                                 "(an object: lattice, atoms, optional name / pp)");
    if (j.at("structure").is_object()) { r.inlineStructure=f.At("structure"); r.structure.clear(); }
    else                               { f.Get("structure",r.structure); r.inlineStructure=json(); }
    if (j.contains("kmesh"))
    {
        const json& k=f.At("kmesh");
        if (!k.is_array() || k.size()!=3) throw std::runtime_error("deck: 'kmesh' must be [nx,ny,nz], got "+k.dump());
        r.kmesh=ivec3_t(k.at(0).get<int>(),k.at(1).get<int>(),k.at(2).get<int>());
    }
    else f.At("kmesh");
    if (j.contains("basis"))
    {
        BasisFromJson(f.At("basis"),r.basis);
    }
    else f.At("basis");
    if (j.contains("state"))
    {
        Fields sf(f.At("state"),"state"); sf.Get("save",r.state.save); sf.Get("restartFrom",r.state.restartFrom); sf.Done();
    }
    else f.At("state");
    if (j.contains("solid")) FromJson(f.At("solid"),r.solid); else f.At("solid");
    if (j.contains("threads"))
    {
        Fields tf(f.At("threads"),"threads"); tf.Get("openmp",r.threads.openmp); tf.Get("blas",r.threads.blas); tf.Done();
        if (r.threads.openmp<0 || r.threads.blas<0)
            throw std::runtime_error("deck: 'threads' counts must be >= 0 (0 = the process's own: QCHEM_OPENMP_THREADS / QCHEM_BLAS_THREADS)");
    }
    else f.At("threads");
    // THE ONE ACCURACY TARGET eps is the RUN's (V-FACADE-ACCURACY, user ruling 2026-11-03): the top-level `accuracy` key, copied into the
    // options of the facade the deck builds (today always the solid's).  `solid.accuracy` is its DEPRECATED alias: read, and warned
    // about.  Both stated is refused -- two statements of one number, and the deck would have to pick.
    const bool solidStates = j.contains("solid") && j.at("solid").is_object() && j.at("solid").contains("accuracy");
    if (j.contains("accuracy"))
    {
        if (solidStates)
            throw std::runtime_error("deck: 'accuracy' and 'solid.accuracy' are both stated -- 'solid.accuracy' is the deprecated spelling of the "
                                     "run's ONE accuracy target; keep the top-level 'accuracy' only");
        double eps=r.solid.accuracy.eps;
        f.Get("accuracy",eps);
        if (!(eps>=0.0)) throw std::runtime_error("deck: 'accuracy' (Ha per atom on an energy run, per matrix element on a geometry run) "
                                                  "must be >= 0 (0 = no target: the rasters tolerances.rasterX ... state)");
        r.solid.accuracy.eps=eps;
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
    // THE XC CHOICES (user ruling 2026-11-09, V-XC-ONE-TERM): `xc.fit` = the Vxc fit basis (was `solid.vxcFit`, refused by name above),
    // `xc.rho` = "exact" (rho[D], the default) | "fit" (the fit's projection P rho).  The run's, copied into the facade's options.
    if (j.contains("xc"))
    {
        Fields xf(f.At("xc"),"xc"); xf.GetEnum("fit",r.solid.vxcFit,kVxc); xf.GetEnum("rho",r.solid.xcRho,kRho); xf.Done();
    }
    else f.At("xc");
    if (j.contains("acceleratorOptions")) FromJson(f.At("acceleratorOptions"),r.accel); else f.At("acceleratorOptions");
    if (j.contains("postSCF"))
    {
        r.postSCF.clear(); size_t i=0;
        for (const auto& a : f.At("postSCF")) r.postSCF.push_back(PostFromJson(a,"postSCF."+std::to_string(i++)));
    }
    else f.At("postSCF");
    if (j.contains("scf") && j.contains("schedule"))
        throw std::runtime_error("deck: give 'scf' (one stage) OR 'schedule' (an annealed recipe), not both -- which would run?");
    if (j.contains("scf")) FromJson(f.At("scf"),r.scf); else f.At("scf");
    if (j.contains("schedule"))
    {
        r.schedule.clear(); size_t i=0;
        for (const auto& st : f.At("schedule"))
        {
            const std::string path="schedule."+std::to_string(i++);
            Fields sf(st,path); RunSpec::Stage stage;
            sf.GetEnum("accelerator",stage.accelerator,kAcc);
            if (st.contains("scf")) FromJson(sf.At("scf"),stage.scf); else sf.At("scf");
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
        try { m=Materials::Inline(spec.inlineStructure.dump(), spec.solid.pp); }
        catch (const std::runtime_error& e) { throw std::runtime_error("deck: structure (inline) -- "+std::string(e.what())); }
    }
    else
    {
        if (StructureData::KindOf(spec.structure)==StructureData::Kind::Molecule)      // KindOf THROWS, listing every known name, on a miss
            throw std::runtime_error("deck: structure '"+spec.structure+"' is a molecule (molecules.json); the deck currently drives "
                                     "periodic runs only (a cell from materials.json)");
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
    if (spec.basis.autoTrim && spec.solid.orthoTol<=0.0)
        throw std::runtime_error("deck: basis.autoTrim needs a solid.orthoTol > 0 to trim at -- at 0 it would drop nothing");
    if (spec.basis.angular==BasisSet::Gaussian::Angular::Spherical && BasisSet::Gaussian::BasisFileOfName(spec.basis.data)=="valence_lowq_sr.bsd")
    {
        bool tm=false; m.cell->ForEachSite([&](int z, const rvec3_t&, bool){ tm = tm || (z>=21 && z<=30); });
        if (tm) throw std::runtime_error("deck: basis.angular \"spherical\" with VALENCE_LOWQ_SR -- the SR transition-metal block's s span lives in the Cartesian d "
                                         "contaminants, which the spherical view removes; use VALENCE_LOWQ_VA or VALENCE_LOWQ_SPH");
    }
    // THE MAGNETIC ORDERING IS STATED (OpenWork §3, 2026-10-26): the facade's own check, run here so a deck is refused BEFORE a revision is claimed.
    try { RequireStatedMagneticOrdering(spec.solid, *m.cell); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: solid.siteSpins -- "+std::string(e.what()).substr(std::string("SolidCalcOptions: ").size())); }
    try { RequireConsistentSiteCharges(spec.solid, *m.cell); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: solid.siteCharges -- "+std::string(e.what()).substr(std::string("SolidCalcOptions: ").size())); }
    // AUTO PICKS THE KIND, NOT THE RECIPE (OpenWork §3, 2026-10-26): the mesh policy's own refusal, run here so --check catches it too.
    try { (void)qcMesh::ResolveXCMesh(spec.solid.xcMesh); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: solid."+std::string(e.what())); }

    // PRE-FLIGHT of the postSCF list (each message names the action and what to change): the old probes found these out AFTER the SCF.
    {
        size_t i=0;
        bool displaced=false;
        for (const auto& a : spec.postSCF)
        {
            const std::string who="deck: postSCF."+std::to_string(i++)+" ("+PostName(a.kind)+"): ";
            if (a.kind==PostAction::Kind::GridAudit)   // the one action that is not about +U: it needs a raster, and a different one
            {
                if (spec.basis.family!=SolidCalcOptions::OrbitalBasis::GPW) throw std::runtime_error(who+"a plane-wave run has no raster to audit");
                if (!(a.referenceMargin>0.0)) throw std::runtime_error(who+"referenceMargin must be > 0 (the reference is FINER than the run)");
                continue;
            }
            if (a.kind==PostAction::Kind::XCMeshAudit) // the run's converged D on a reference Becke mesh: a Gaussian run with a positive mesh
            {
                if (spec.basis.family!=SolidCalcOptions::OrbitalBasis::GPW) throw std::runtime_error(who+"a plane-wave run samples XC on its FFT grid, not a Becke mesh");
                if (a.xcRefNRadial<=0 || a.xcRefAngularDegree<=0) throw std::runtime_error(who+"nRadial and angularDegree must be > 0");
                // displace runs in its own session, so this run's D stays on the deck's positions: after one, the audit would read as the displaced point's.
                if (displaced) throw std::runtime_error(who+"it audits the deck's own positions, not a displaced point -- list it before displace, or state the displaced structure in the deck");
                continue;
            }
            if (a.kind==PostAction::Kind::EggBox)      // a rigid translation of a free, U = 0 Gaussian run at a nonzero shift
            {
                if (spec.basis.family!=SolidCalcOptions::OrbitalBasis::GPW) throw std::runtime_error(who+"a plane-wave basis is translation-invariant by construction");
                if (spec.solid.imposeSymmetry) throw std::runtime_error(who+"an imposed group's SALC phases follow the origin -- gate the free twin (imposeSymmetry false)");
                if (!spec.solid.hubbard.empty()) throw std::runtime_error(who+"the +U projectors would not follow the shift -- gate the U = 0 twin");
                if (norm(a.shift)==0.0) throw std::runtime_error(who+"a zero shift measures nothing");
                continue;
            }
            if (a.kind==PostAction::Kind::Displace)    // geometry points through ONE SolidSession (its own setup, sized on max|dH_ij|)
            {
                if (a.points.empty()) throw std::runtime_error(who+"no points -- list {\"atoms\":[i,...], \"delta\":[x,y,z]} (bohr)");
                displaced=true;
                // The session's own refusal, run here so a deck fails before its SCF: a free mesh made invariant under the group detected at
                // EACH geometry switches along a path that breaks it (OpenWork, closed 2026-10-30).
                if (!spec.solid.imposeSymmetry && spec.basis.family==SolidCalcOptions::OrbitalBasis::GPW && spec.solid.xcMeshSymmetric)
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
            if (needsU && spec.solid.hubbard.empty())
                throw std::runtime_error(who+"the run carries no Hubbard manifold -- list the channels in solid.hubbard (at U_eV=0 to probe without +U)");
            const bool lrLoop = a.kind==PostAction::Kind::HubbardLoop && a.estimator==SolidCalculation::HubbardLoop::Estimator::LinearResponse;
            const bool isResponse = a.kind==PostAction::Kind::IndependentResponse || a.kind==PostAction::Kind::HubbardLinearResponse || lrLoop;
            // NOT an error: whether an imposed group REDUCES the k-mesh depends on the group and the mesh (the NiO AFM-II 2x2x2 Shubnikov group maps 8 k-points onto 8
            // irreducible ones, and its chi0 ran).  A blanket refusal would reject decks the library accepts; the response itself THROWS at run time if the mesh was
            // really reduced (D5).  Vet() says so in a note.
            // A smeared (metal) response at q != 0 needs ONE Fermi level across the k-mesh (Reference::Gap refuses it at run time --
            // after the SCF; SrVO3 2026-10-05 ran two hours to find that out).  The FINAL stage's kT is the response's occupancy rule.
            if (isResponse && spec.Stages().back().scf.SmearingkT>0.0 && spec.solid.perKFermiLevels)
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
            if ((a.kind==PostAction::Kind::HubbardLinearResponse || lrLoop) && !spec.solid.forceComplex)
                throw std::runtime_error(who+"needs the complex ansatz (the real-TRIM response face is not built); set solid.forceComplex true");
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
                if (j>=spec.solid.hubbard.size())
                    throw std::runtime_error(who+"perturb index "+std::to_string(j)+" is past the "+std::to_string(spec.solid.hubbard.size())+" manifolds in solid.hubbard");
        }
    }
    // The material NAMES its pseudopotential set; a deck that left `solid.pp` at the facade default takes the material's.  `species` (the
    // per-element overrides) is filled with the resolved vocabulary so the record states every valence actually used.
    if (spec.solid.pp==SolidCalcOptions{}.pp) spec.solid.pp=m.pp;
    if (spec.solid.species.empty()) spec.solid.species=m.species;
    // The electron count of the species the RUN uses: a stated override (Na q9) changes it.  Until 2026-10-31 this read the
    // material's own species, so a deck stating Na q9 ran NaF with 8 electrons instead of 16 -- a cell charged +8, silently.
    if (spec.solid.Nelec==0)
    {
        Materials::Material used=m;
        for (auto& [el,q] : used.species)
            for (const auto& [el2,q2] : spec.solid.species) if (el2==el) q=q2;
        spec.solid.Nelec=used.Nelec()-spec.solid.charge;   // a STATED net charge (SolidCalcOptions::charge)
    }
    // THE SPIN BOOKKEEPING on the resolved count (the setup's own decision 2, run here so --check refuses an impossible request:
    // multiplicity 1 on an odd count once ran a silent doublet).
    try { (void)SpinChannelCounts(spec.solid, spec.solid.Nelec); }
    catch (const std::invalid_argument& e) { throw std::runtime_error("deck: solid.multiplicity -- "+std::string(e.what()).substr(std::string("SolidCalcOptions: ").size())); }
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
            try { (void)qchem::ResolveMixing(stages[i].scf, /*gSpace*/true, fock); }
            catch (const std::invalid_argument& e)
            {
                const std::string where = spec.schedule.empty() ? std::string("scf") : "schedule."+std::to_string(i)+".scf";
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


VetReport Vet(RunSpec spec, const std::filesystem::path& outDir)
{
    namespace fs=std::filesystem;
    const Materials::Material mat=Resolve(spec);
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
    std::shared_ptr<const BasisSet::Real_BS> b=recipe.Make(*mat.cell, Pseudopotential::PseudopotentialSet::Get(spec.solid.pp).With(spec.solid.species));
    size_t nf=0;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*b).Iterate<BasisSet::Real_OIBS>()) nf+=ibs->GetNumFunctions();
    {
        bool resp=false; for (const auto& a : spec.postSCF) resp = resp || a.kind==RunSpec::PostAction::Kind::IndependentResponse || a.kind==RunSpec::PostAction::Kind::HubbardLinearResponse;
        if (resp && spec.solid.imposeSymmetry && spec.kmesh.x*spec.kmesh.y*spec.kmesh.z>1)
            v.notes.push_back("a response with solid.imposeSymmetry on a k-mesh: if the imposed group REDUCES the mesh the response throws at run time (D5); the NiO AFM-II 2x2x2 Shubnikov group does not reduce it");
        if (resp && LooserThan(spec.Stages().back().scf.MinΔρPerElectron, kResponseReferenceMinΔρ))
            v.notes.push_back("a response reference should converge to scf.minDeltaRhoPerElectron <= 1e-9 (a few more iterations; ρ_lost/N does not floor it): "
                              "at 1e-6 its orbitals carry ~1e-6 symmetry and eigenpair residuals");
    }
    if (spec.basis.autoTrim) v.notes.push_back("basis.autoTrim: the vet-stage trim loop is NOT run by --check (function count is the untrimmed one)");
    std::ostringstream os;
    os<<spec.StructureName()<<(spec.inlineStructure.is_null() ? "" : " (inline)")<<": "<<mat.cell->GetNumAtoms()<<" atoms, Nelec "<<spec.solid.Nelec<<", multiplicity "<<spec.solid.multiplicity
      <<", k-mesh "<<spec.kmesh.x<<"x"<<spec.kmesh.y<<"x"<<spec.kmesh.z<<", basis "<<(!spec.basis.inlineBsd.empty() ? spec.basis.data+" (inline)" : spec.basis.file.empty()?spec.basis.data:spec.basis.file)<<(spec.basis.angular==BasisSet::Gaussian::Angular::Spherical?" (spherical)":"")
      <<" "<<nf<<" functions, "<<spec.Stages().size()<<" SCF stage(s), "<<spec.solid.hubbard.size()<<" Hubbard manifold(s), "<<spec.postSCF.size()<<" postSCF action(s)"
      <<(spec.solid.imposeSymmetry?", symmetry IMPOSED":", symmetry FREE");
    v.summary=os.str();
    return v;
}

// postSCF displace: the deck's own setup (re-resolved: an UNSIZED one, so the session's reference pilot sizes it on max|dH_ij|), the
// reference at the deck's positions, then each point -- all through ONE SolidSession, so every point shares the reference's frozen
// basis recipe and rasters.  The main run is not reused: it was sized as an Energy.  No saved state (the main run's file is the deck's).
static std::shared_ptr<const SolidSetup> SetupOf(const RunSpec&, const Materials::Material&, const std::vector<SCFStage>&);
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
    spec.solid.orbitalBasis=spec.basis.family; spec.solid.planeWaveEcut=spec.basis.ecut;
    spec.solid.saveStateTo.clear();
    std::vector<SCFStage> stages;
    for (const auto& st : spec.Stages()) stages.push_back({st.scf,st.accelerator});   // the whole recipe: the session seeds, it never restarts
    SolidSession S(SetupOf(spec, mat, stages), stages);
    const PESGeometry ref=S.Reference();
    DisplaceResult out;
    std::cout<<"[displace] reference: the deck's positions ("<<ref.R.size()<<" atoms), one SolidSession for "<<moves.size()<<" point(s)"<<std::endl;
    auto p0=S.ConvergeAt(ref.cell, ref.R);
    if (!p0) { out.failure0=p0.Error().details; return out; }
    out.E0=p0->energy; out.it0=p0->iterations;
    bool ok=true;
    for (size_t k=0;k<moves.size();++k)
    {
        std::vector<rvec3_t> R=ref.R;
        if (moves[k].atoms.empty()) for (auto& r : R) r+=moves[k].delta;
        else                        for (size_t a : moves[k].atoms) R[a]+=moves[k].delta;
        auto p=S.ConvergeAt(ref.cell, R);
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
    prov.policyResolved=RunPolicy(spec.solid.policy).Banner();   // the stated policy is in the deck; what it RESOLVED to is the record
    prov.basisFile=!spec.basis.inlineBsd.empty() ? MaterializedBasis(spec.basis.inlineBsd)
                  : spec.basis.file.empty() ? BasisSet::Gaussian::BasisFileOfName(spec.basis.data) : spec.basis.file;
    prov.basisChecksum=BasisSet::Gaussian::BasisFileChecksum(prov.basisFile);
    WriteRevision(out.revision, ToJson(spec), prov);     // the record exists BEFORE the SCF: a crashed run still leaves it
    const std::string stem=out.revision.stem().string(); // <structure>.rNNN -- the run's name in its own output

    if (spec.solid.label=="gpw") spec.solid.label=stem;
    if (!spec.state.save.empty())
    {
        if (spec.state.save=="auto") { fs::create_directories(outDir/"states"); spec.state.save=(outDir/"states"/(stem+".h5")).string(); }
    }
    if (!restartPath.empty()) spec.state.restartFrom=restartPath;   // Build takes a path: the stem was this outDir's
    const std::vector<SCFStage> stages=ScheduleToRun(spec);

    qchem::report::Begin(stem);
    qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
    const Close close;
    const std::unique_ptr<SolidCalculation> calc=Build(spec);   // ONE path with the GUI's: built, then converged here
    (void)calc->Converge(stages);
    if (auto R=calc->Result())
    {
        out.converged=true; out.energy=calc->LastIterateTerms().GetTotalEnergy();
        out.summary="CONVERGED  "+calc->Diagnostics().Summary();
    }
    else out.summary=(R.Error().why==SCFFailure::Why::ResourceRefused ? "REFUSED ON RESOURCES: " : "NOT converged: ")+R.Error().details;
    // A run REFUSED ON RESOURCES (the audit RSS policy) did not reach the state its postSCF actions read: before the first SCF there is no D
    // at all, and mid-pilot the rasters are not settled.  Each action is recorded as skipped, and the record carries the refusal's numbers.
    std::optional<ResourceRefusal> refused;
    if (const auto R=calc->Result(); !R && R.Error().resource) refused=R.Error().resource;

    // postSCF: each action on the converged calculation, in order, each reporting itself; the summaries go in the revision's `results`.  A failed
    // action (a response gate, an unconverged FD step) is a RESULT, not an exception: the remaining actions still run.
    json results=json::array();
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
                return qcMesh::XCAuditReference(S->xcMesh, a.xcRefNRadial, a.xcRefAngularDegree, S->opts.beckeLadder);
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
            const double target=S->opts.accuracy.Target(c, S->cell->GetNumAtoms());
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
    qchem::report::EmitTimings();   // sorted by cost + PEAK RSS, inside the report bracket (what scripts/bench greps)
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
    for (const auto& s : spec.Stages()) stages.push_back({s.scf,s.accelerator});
    if (!spec.state.restartFrom.empty()) stages={stages.back()};   // a restart is ONE stage: the schedule's final one
    return stages;
}

// THE SETUP, resolved ONCE: the basis from the deck's recipe (the trim stated, or vetted at the run's own orthoTol and frozen), and every
// decision the run makes before it builds anything.  The first stage's accelerator (a restart's: its one stage's) is the setup's; each
// later stage names its own.  \a spec is Build's (its family copied into solid, its save path resolved).
static std::shared_ptr<const SolidSetup> SetupOf(const RunSpec& spec, const Materials::Material& mat, const std::vector<SCFStage>& stages)
{
    Lattice_3D lat(*mat.cell, spec.kmesh);
    SolidCalcOptions so=spec.solid; so.accelerator=stages.front().accelerator;
    return SolidSetup::Resolve(lat, RecipeOf(spec.basis), so, spec.accel);
}

std::unique_ptr<SolidCalculation> Build(const RunSpec& specIn, SolidCalculation::Observer onIteration)
{
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
    spec.solid.orbitalBasis=spec.basis.family; spec.solid.planeWaveEcut=spec.basis.ecut;   // the basis block names the family; the facade reads it from solid
    if (!spec.state.save.empty())
    {
        if (spec.state.save=="auto") throw std::runtime_error("deck: state.save 'auto' names a file under Run's outDir; Build takes a path");
        spec.solid.saveStateTo=spec.state.save;
    }
    const std::vector<SCFStage> stages=ScheduleToRun(spec);
    const auto setup=SetupOf(spec, mat, stages);
    if (r.empty()) return SolidCalculation::Prepare(setup, stages, std::move(onIteration));
    // a refusal THROWS -- a deck never silently re-seeds
    auto R=SolidCalculation::PrepareRestart(r, setup, stages.front().params, std::move(onIteration));
    if (!R) throw std::runtime_error("deck: state.restartFrom "+r+" refused: "+R.Error().details);
    return R.TakeValue();
}

//=== The shell what-if: the deck's setup, its UNTRIMMED basis, and the shells its trim keeps =======================================
BasisShells MakeBasisShells(RunSpec spec)
{
    const Materials::Material mat=Resolve(spec);
    spec.solid.orbitalBasis=spec.basis.family; spec.solid.planeWaveEcut=spec.basis.ecut;   // as Build: the basis block names the family
    BasisShells out;
    out.setup=SetupOf(spec, mat, ScheduleToRun(spec));
    out.table=out.setup->MakeShellConditioning();
    out.floor=spec.solid.orthoTol>0.0 ? spec.solid.orthoTol : ConditioningFloor;
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
