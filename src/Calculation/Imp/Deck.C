// File: Calculation/Imp/Deck.C  See the interface.
module;
#include <cstdio>
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

//=== SolidCalcOptions ===========================================================================
namespace
{
json ToJson(const Hamiltonian::HubbardManifold& h)
{
    return {{"site",h.site},{"l",h.l},{"U_Ha",h.U},{"UirrepHa",h.Uirrep},{"alphaHa",h.alpha},{"radial",h.radial},
            {"atomicRadial",h.atomicRadial},{"orthoAtomic",h.orthoAtomic}};
}
void FromJson(const json& j, Hamiltonian::HubbardManifold& h, const std::string& path)
{
    Fields f(j,path);
    f.Get("site",h.site); f.Get("l",h.l); f.Get("U_Ha",h.U); f.Get("UirrepHa",h.Uirrep); f.Get("alphaHa",h.alpha);
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
            {"xcMesh",ToJson(o.xcMesh)},{"tolerances",ToJson(o.tolerances)},{"vxcFit",NameOf(o.vxcFit,kVxc)},{"hubbard",hub},
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
json ToJson(const RunSpec& r) { return {{"structure",r.structure},{"solid",ToJson(r.solid)},{"scf",ToJson(r.scf)}}; }
void FromJson(const json& j, RunSpec& r)
{
    Fields f(j,"run");
    if (!j.contains("structure") || !j.at("structure").is_string())
        throw std::runtime_error("deck: 'structure' is required and must be a NAME from materials.json or molecules.json");
    f.Get("structure",r.structure);
    if (j.contains("solid")) FromJson(f.At("solid"),r.solid); else f.At("solid");
    if (j.contains("scf"))   FromJson(f.At("scf"),r.scf);     else f.At("scf");
    f.Done();
}
Materials::Material Resolve(RunSpec& spec)
{
    if (StructureData::KindOf(spec.structure)==StructureData::Kind::Molecule)      // KindOf THROWS, listing every known name, on a miss
        throw std::runtime_error("deck: structure '"+spec.structure+"' is a molecule (molecules.json); the deck currently drives "
                                 "periodic runs only (a cell from materials.json)");
    Materials::Material m=Materials::Get(spec.structure);
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
    json pv={{"commandLine",prov.commandLine},{"overrides",prov.overrides},{"ignoredEnvironment",prov.ignoredEnvironment}};
    if (!prov.inputDeck.empty())
    {
        pv["inputDeck"]=prov.inputDeck.string();
        if (std::filesystem::exists(prov.inputDeck)) pv["inputDeckChecksum"]=Fnv1a(prov.inputDeck);
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

} // namespace qchem::deck
