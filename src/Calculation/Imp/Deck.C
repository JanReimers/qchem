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
#include <fstream>
#include <functional>
#include <initializer_list>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>
module qchem.Deck;

import qchem.BasisSet.Lattice.BasisSet;    // RasterPolicy, CellImages
import qchem.Hamiltonian.Factory;           // VxcFit, HubbardManifold
import qchem.SCFAccelerator.Factory;        // SCFAccelerators::Type
import qchem.ChargeDensity.Seed;            // SeedStrategy
import qchem.LASolver;                      // qchem::Ortho
import qchem.Types;
import qchem.Materials;
import qchem.SolidCalculation;
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;
import qchem.Reporting;
import qchem.Outcome;
import qchem.BasisSet.Lattice.BasisSet;   // VetStageTrim, GPWParams
import qchem.BasisSet.Gaussian.Point.ShellTrim;
import qchem.RunPolicy;

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
            if (!itsKnown.count(it.key()))
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

using BasisSet::PlaneWave::RasterPolicy;
using BasisSet::Gaussian::CellImages;
using ChargeDensity::SeedStrategy;
using AccType = SCFAccelerators::Type;
using qcMesh::UnitCellKind; using qcMesh::RadialKind; using qcMesh::AngularKind;

#define NAMES(T,...) std::initializer_list<std::pair<const char*,T>>{__VA_ARGS__}
const auto kRaster = NAMES(RasterPolicy, {"AliasFree",RasterPolicy::AliasFree}, {"BallOnly",RasterPolicy::BallOnly});
const auto kImages = NAMES(CellImages,   {"Periodic",CellImages::Periodic}, {"HomeCellOnly",CellImages::HomeCellOnly});
const auto kVxc    = NAMES(Hamiltonian::VxcFit, {"Auto",Hamiltonian::VxcFit::Auto}, {"PlaneWave",Hamiltonian::VxcFit::PlaneWave}, {"Delta",Hamiltonian::VxcFit::Delta});
const auto kAcc    = NAMES(AccType, {"DIIS",AccType::DIIS}, {"GDM",AccType::GDM}, {"Ladder",AccType::Ladder}, {"Null",AccType::Null});
const auto kSeed   = NAMES(SeedStrategy, {"Default",SeedStrategy::Default}, {"CoreGuess",SeedStrategy::CoreGuess}, {"Uniform",SeedStrategy::Uniform},
                                         {"SAD",SeedStrategy::SAD}, {"IonicSAD",SeedStrategy::IonicSAD});
const auto kOrtho  = NAMES(Ortho, {"Cholesky",Cholesky}, {"Eigen",Eigen}, {"SVD",SVD}, {"Auto",Auto}, {"CholeskyPivoted",CholeskyPivoted});
const auto kCell   = NAMES(UnitCellKind, {"Uniform",UnitCellKind::Uniform}, {"Becke",UnitCellKind::Becke}, {"Auto",UnitCellKind::Auto});
const auto kRadial = NAMES(RadialKind, {"MHL",RadialKind::MHL}, {"Log",RadialKind::Log}, {"Linear",RadialKind::Linear});
const auto kAngular= NAMES(AngularKind, {"Lebedev",AngularKind::Lebedev}, {"GaussLegendre",AngularKind::GaussLegendre});
using Measure = SCFParams::Measure;
using BasisSet::Gaussian::BasisSetData;
const auto kBasisData = NAMES(BasisSetData, {"DZVP",BasisSetData::DZVP}, {"DZVP2",BasisSetData::DZVP2}, {"TZVP",BasisSetData::TZVP}, {"ORB",BasisSetData::ORB},
    {"ORB1",BasisSetData::ORB1}, {"SIPP",BasisSetData::SIPP}, {"SIPP_SR",BasisSetData::SIPP_SR}, {"VALENCE_LOWQ",BasisSetData::VALENCE_LOWQ},
    {"VALENCE_LOWQ_SR",BasisSetData::VALENCE_LOWQ_SR}, {"VALENCE_LOWQ_SR2",BasisSetData::VALENCE_LOWQ_SR2}, {"VALENCE_LOWQ_SPH",BasisSetData::VALENCE_LOWQ_SPH},
    {"VALENCE_LOWQ_VA",BasisSetData::VALENCE_LOWQ_VA}, {"VALENCE_LOWQ_VB",BasisSetData::VALENCE_LOWQ_VB});
const auto kMeasure= NAMES(Measure, {"MixerResidual",Measure::MixerResidual}, {"MaxDeltaD",Measure::MaxΔD});

} // anon

//=== GPWTolerances ==============================================================================
json ToJson(const BasisSet::Gaussian::GPWTolerances& t)
{
    return {{"vlocEps",t.vlocEps},{"localPPRelCutoff",t.localPPRelCutoff},{"relFieldSharp",t.relFieldSharp},{"mgridEcuts",t.mgridEcuts},
            {"screenEps",t.screenEps},{"fieldSharp",t.fieldSharp},{"relCutoff",t.relCutoff},{"densityEps",t.densityEps}};
}
void FromJson(const json& j, BasisSet::Gaussian::GPWTolerances& t)
{
    Fields f(j,"tolerances");
    f.Get("vlocEps",t.vlocEps); f.Get("localPPRelCutoff",t.localPPRelCutoff); f.Get("relFieldSharp",t.relFieldSharp);
    f.Get("mgridEcuts",t.mgridEcuts); f.Get("screenEps",t.screenEps); f.Get("fieldSharp",t.fieldSharp);
    f.Get("relCutoff",t.relCutoff); f.Get("densityEps",t.densityEps);
    f.Done();
}

//=== SCFParams ==================================================================================
json ToJson(const SCFParams& p)
{
    return {{"NMaxIter",p.NMaxIter},{"minDeltaRho",p.MinΔρ},{"deltaRhoMeasure",NameOf(p.Δρmeasure,kMeasure)},{"minDeltaFD",p.MinΔFD},
            {"minDeltaE",p.MinΔE},{"minVirial",p.MinVirial},{"minFD",p.MinFD},{"startingRelaxRo",p.StartingRelaxRo},{"mergeTol",p.MergeTol},
            {"verbose",p.Verbose},{"xcCuspDeficit",p.XCCuspDeficit},{"kerkerG0",p.KerkerG0},{"useMOM",p.UseMOM},{"pulayDepth",p.PulayDepth},
            {"pulayStart",p.PulayStart},{"momStartIter",p.MOMStartIter},{"smearingkT",p.SmearingkT},{"momSmearPenalty",p.MOMSmearPenalty},
            {"stopOnAccelExhausted",p.StopOnAccelExhausted},
            {"momGuard",{{"holePersistence",p.Guard.HolePersistence},{"maxReleases",p.Guard.MaxReleases}}}};
}
void FromJson(const json& j, SCFParams& p)
{
    Fields f(j,"scf");
    f.Get("NMaxIter",p.NMaxIter); f.Get("minDeltaRho",p.MinΔρ); f.GetEnum("deltaRhoMeasure",p.Δρmeasure,kMeasure);
    f.Get("minDeltaFD",p.MinΔFD); f.Get("minDeltaE",p.MinΔE); f.Get("minVirial",p.MinVirial); f.Get("minFD",p.MinFD);
    f.Get("startingRelaxRo",p.StartingRelaxRo); f.Get("mergeTol",p.MergeTol); f.Get("verbose",p.Verbose);
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

//=== MeshParams =================================================================================
json ToJson(const qcMesh::MeshParams& m)
{
    return {{"radial",NameOf(m.radial,kRadial)},{"nRadial",m.nRadial},{"mhlM",m.mhl_m},{"mhlAlpha",m.mhl_alpha},{"logStart",m.logStart},
            {"logStop",m.logStop},{"angular",NameOf(m.angular,kAngular)},{"angularDegree",m.angularDegree},{"angRot",m.angRot},
            {"beckeOrder",m.beckeOrder},{"nUniform",m.nUniform},{"eCut",m.eCut},{"relCutoff",m.relCutoff},
            {"cellKind",NameOf(m.cellKind,kCell)},{"beckeEps",m.beckeEps}};
}
void FromJson(const json& j, qcMesh::MeshParams& m)
{
    Fields f(j,"xcMesh");
    f.GetEnum("radial",m.radial,kRadial); f.Get("nRadial",m.nRadial); f.Get("mhlM",m.mhl_m); f.Get("mhlAlpha",m.mhl_alpha);
    f.Get("logStart",m.logStart); f.Get("logStop",m.logStop); f.GetEnum("angular",m.angular,kAngular);
    f.Get("angularDegree",m.angularDegree); f.Get("angRot",m.angRot); f.Get("beckeOrder",m.beckeOrder); f.Get("nUniform",m.nUniform);
    f.Get("eCut",m.eCut); f.Get("relCutoff",m.relCutoff); f.GetEnum("cellKind",m.cellKind,kCell); f.Get("beckeEps",m.beckeEps);
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
    return {{"site",h.site},{"l",h.l},{"U_eV",ToEV(h.U)},{"Uirrep_eV",ToEV(h.Uirrep)},{"alpha_eV",ToEV(h.alpha)},{"radial",h.radial},
            {"atomicRadial",h.atomicRadial},{"orthoAtomic",h.orthoAtomic}};
}
void FromJson(const json& j, Hamiltonian::HubbardManifold& h, const std::string& path)
{
    Fields f(j,path);
    f.Get("site",h.site); f.Get("l",h.l);
    double u=h.U*kHaToEV; f.Get("U_eV",u); h.U=u/kHaToEV;                                   // eV in the file, a.u. in RAM
    std::vector<double> ui=ToEV(h.Uirrep); f.Get("Uirrep_eV",ui); h.Uirrep=FromEV(ui);
    double al=h.alpha*kHaToEV; f.Get("alpha_eV",al); h.alpha=al/kHaToEV;
    f.Get("radial",h.radial); f.Get("atomicRadial",h.atomicRadial); f.Get("orthoAtomic",h.orthoAtomic);
    f.Done();
}
} // anon

json ToJson(const SolidCalcOptions& o)
{
    json species=json::array();
    for (const auto& s : o.species) species.push_back({s.first,s.second});
    json hub=json::array();
    for (const auto& h : o.hubbard) hub.push_back(ToJson(h));
    return {{"Nelec",o.Nelec},{"multiplicity",o.multiplicity},{"species",species},
            {"densityEcut",o.densityEcut},{"cutoffFactor",o.cutoffFactor},{"ladderFactor",o.ladderFactor},
            {"raster",NameOf(o.raster,kRaster)},{"images",NameOf(o.images,kImages)},{"kShift",{o.kShift.x,o.kShift.y,o.kShift.z}},
            {"xcMesh",ToJson(o.xcMesh)},{"tolerances",ToJson(o.tolerances)},{"policy",ToJson(o.policy)},{"vxcFit",NameOf(o.vxcFit,kVxc)},{"hubbard",hub},
            {"accelerator",NameOf(o.accelerator,kAcc)},{"globalFermi",o.globalFermi},{"imposeSymmetry",o.imposeSymmetry},
            {"seed",NameOf(o.seed,kSeed)},{"ortho",NameOf(o.ortho,kOrtho)},{"orthoTol",o.orthoTol},{"forceComplex",o.forceComplex},
            {"spinsShareFermi",o.spinsShareFermi},{"greyImposition",o.greyImposition},{"momFromSeed",o.momFromSeed},
            {"siteSpins",o.siteSpins},{"label",o.label},{"saveStateTo",o.saveStateTo}};
}
void FromJson(const json& j, SolidCalcOptions& o)
{
    Fields f(j,"solid");
    f.Get("Nelec",o.Nelec); f.Get("multiplicity",o.multiplicity);
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
    f.GetEnum("raster",o.raster,kRaster); f.GetEnum("images",o.images,kImages);
    if (j.contains("kShift"))
    {
        const json& k=f.At("kShift");
        if (!k.is_array() || k.size()!=3) throw std::runtime_error("deck: 'solid.kShift' must be [x,y,z], got "+k.dump());
        o.kShift=rvec3_t(k.at(0).get<double>(),k.at(1).get<double>(),k.at(2).get<double>());
    }
    else f.At("kShift");
    if (j.contains("xcMesh"))     FromJson(f.At("xcMesh"),o.xcMesh);         else f.At("xcMesh");
    if (j.contains("tolerances")) FromJson(f.At("tolerances"),o.tolerances); else f.At("tolerances");
    if (j.contains("policy"))     FromJson(f.At("policy"),o.policy);         else f.At("policy");
    f.GetEnum("vxcFit",o.vxcFit,kVxc);
    if (j.contains("hubbard"))
    {
        o.hubbard.clear();
        size_t i=0;
        for (const auto& h : f.At("hubbard")) { o.hubbard.emplace_back(); FromJson(h,o.hubbard.back(),"solid.hubbard."+std::to_string(i++)); }
    }
    else f.At("hubbard");
    f.GetEnum("accelerator",o.accelerator,kAcc); f.Get("globalFermi",o.globalFermi); f.Get("imposeSymmetry",o.imposeSymmetry);
    f.GetEnum("seed",o.seed,kSeed); f.GetEnum("ortho",o.ortho,kOrtho); f.Get("orthoTol",o.orthoTol); f.Get("forceComplex",o.forceComplex);
    f.Get("spinsShareFermi",o.spinsShareFermi); f.Get("greyImposition",o.greyImposition); f.Get("momFromSeed",o.momFromSeed);
    f.Get("siteSpins",o.siteSpins); f.Get("label",o.label); f.Get("saveStateTo",o.saveStateTo);
    f.Done();
}

//=== RunSpec ====================================================================================
namespace
{
json BasisJson(const RunSpec::Basis& b)
{
    json j={{"data",NameOf(b.data,kBasisData)},{"spherical",b.spherical}};
    if (b.vet) j["vet"]=true;
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
    f.GetEnum("data",b.data,kBasisData); f.Get("spherical",b.spherical); f.Get("vet",b.vet);
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
    json j={{"structure",r.structure},{"kmesh",{r.kmesh.x,r.kmesh.y,r.kmesh.z}},
            {"basis",BasisJson(r.basis)},{"solid",ToJson(r.solid)}};
    if (!r.state.save.empty() || !r.state.restartFrom.empty())
    {
        json st=json::object();
        if (!r.state.save.empty())        st["save"]=r.state.save;
        if (!r.state.restartFrom.empty()) st["restartFrom"]=r.state.restartFrom;
        j["state"]=st;
    }
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
    if (!j.contains("structure") || !j.at("structure").is_string())
        throw std::runtime_error("deck: 'structure' is required and must be a NAME from materials.json or molecules.json");
    f.Get("structure",r.structure);
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
    if (StructureData::KindOf(spec.structure)==StructureData::Kind::Molecule)      // KindOf THROWS, listing every known name, on a miss
        throw std::runtime_error("deck: structure '"+spec.structure+"' is a molecule (molecules.json); the deck currently drives "
                                 "periodic runs only (a cell from materials.json)");
    Materials::Material m=Materials::Get(spec.structure);
    if (spec.basis.vet && !spec.basis.trim.empty())
        throw std::runtime_error("deck: basis.vet and basis.trim are exclusive (a vet computes the trim; a stated trim is the A/B for it)");
    if (spec.basis.spherical && spec.basis.data==BasisSetData::VALENCE_LOWQ_SR)
    {
        bool tm=false; m.cell->ForEachSite([&](int z, const rvec3_t&, bool){ tm = tm || (z>=21 && z<=30); });
        if (tm) throw std::runtime_error("deck: basis.spherical with VALENCE_LOWQ_SR -- the SR transition-metal block's s span lives in the Cartesian d "
                                         "contaminants, which the spherical view removes; use VALENCE_LOWQ_VA or VALENCE_LOWQ_SPH");
    }
    if (spec.solid.species.empty()) spec.solid.species=m.species;
    if (spec.solid.Nelec==0)        spec.solid.Nelec=m.Nelec();
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

json LoadDeck(const std::filesystem::path& path, const std::string& currentCodeVersion)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("deck: cannot read "+path.string());
    json j;
    try { j=json::parse(in); } catch (const json::exception& e) { throw std::runtime_error("deck: "+path.string()+": "+e.what()); }
    if (!j.contains("run")) return j;                                  // a hand-written deck: the payload IS the file
    const json& h=j.value("deck",json::object());
    if (h.value("schema",0) > kSchemaVersion)
        throw std::runtime_error("deck: "+path.string()+" has schema "+std::to_string(h.value("schema",0))+", newer than this code's "+std::to_string(kSchemaVersion));
    if (h.value("codeVersion",std::string()) != currentCodeVersion)
        std::fprintf(stderr,"[deck] WARNING: %s was written by code version '%s', this is '%s' -- the run may not reproduce\n",
                     path.c_str(), h.value("codeVersion",std::string()).c_str(), currentCodeVersion.c_str());
    return j.at("run");
}


RunOutcome Run(RunSpec spec, Provenance prov, const std::filesystem::path& outDir)
{
    namespace fs=std::filesystem;
    const Materials::Material mat=Resolve(spec);

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

    RunOutcome out;
    out.revision=ClaimRevision(outDir, spec.structure);
    prov.policyResolved=RunPolicy(spec.solid.policy).Banner();   // the stated policy is in the deck; what it RESOLVED to is the record
    WriteRevision(out.revision, ToJson(spec), prov);     // the record exists BEFORE the SCF: a crashed run still leaves it
    const std::string stem=out.revision.stem().string(); // <structure>.rNNN -- the run's name in its own output

    if (spec.solid.label=="gpw") spec.solid.label=stem;
    if (!spec.state.save.empty())
    {
        if (spec.state.save=="auto") { fs::create_directories(outDir/"states"); spec.solid.saveStateTo=(outDir/"states"/(stem+".h5")).string(); }
        else spec.solid.saveStateTo=spec.state.save;
    }
    Lattice_3D lat(*mat.cell, spec.kmesh);
    // THE BASIS: the data file, optionally trimmed (stated, or by the vet loop at the run's own orthoTol), optionally in the spherical lattice view.
    auto make=[&](const BasisSet::Gaussian::ShellTrim& trim)
    {
        std::shared_ptr<const BasisSet::Real_BS> b(BasisSet::Gaussian::Factory(spec.basis.data, mat.cell.get(),
                                                   BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian, trim));
        if (spec.basis.spherical) b=BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(std::move(b));
        return b;
    };
    std::shared_ptr<const BasisSet::Real_BS> mol;
    if (spec.basis.vet)
        mol=BasisSet::Lattice::VetStageTrim(lat, make, {.images=spec.solid.images, .kShift=spec.solid.kShift}, spec.solid.orthoTol).mol;
    else
    {
        BasisSet::Gaussian::ShellTrim trim;
        for (const auto& t : spec.basis.trim) trim.shells.push_back({t.Z, t.l, rvec_t(1,t.alpha)});
        if (!trim.empty()) { std::cout<<"[basis trim] STATED (no vet loop): "; trim.Write(std::cout); std::cout<<std::endl; }
        mol=make(trim);
    }

    std::vector<SCFStage> stages;
    for (const auto& s : spec.Stages()) stages.push_back({s.scf,s.accelerator});

    qchem::report::Begin(stem);
    qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
    struct Close { ~Close() { qchem::report::ClearConsole(); qchem::report::End(); } } close;
    std::unique_ptr<SolidCalculation> calc;
    if (restartPath.empty()) calc=std::make_unique<SolidCalculation>(lat, mol, spec.solid, stages);
    else
    {   // one stage: the schedule's FINAL one (its params and accelerator); a refusal THROWS -- a deck run never silently re-seeds
        SolidCalcOptions ro=spec.solid; ro.accelerator=stages.back().accelerator;
        auto R=SolidCalculation::Restart(restartPath, lat, mol, ro, stages.back().params);
        if (!R) throw std::runtime_error("deck: state.restartFrom "+restartPath+" refused: "+R.Error().details);
        calc=R.TakeValue();
    }
    if (auto R=calc->Result())
    {
        out.converged=true; out.energy=calc->LastIterateTerms().GetTotalEnergy();
        out.summary="CONVERGED  "+calc->Diagnostics().Summary();
    }
    else out.summary="NOT converged: "+R.Error().details;
    return out;
}

} // namespace qchem::deck
