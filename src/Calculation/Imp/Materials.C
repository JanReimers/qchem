// File: Calculation/Imp/Materials.C  The JSON readers behind qchem.Materials.
module;
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <nlohmann/json.hpp>
module qchem.Materials;
import qchem.PeriodicTable;   // thePeriodicTable().GetZ(symbol)
import qchem.Matrix3D;
import qchem.Types;

namespace qchem::Materials
{

namespace
{
using json = nlohmann::ordered_json;   // ORDERED: file order is the pick-list order

const json& Database(const char* file)
{
    // One parse per file per process; the files are small and read-only.
    static std::map<std::string, json> dbs;
    auto it=dbs.find(file);
    if (it!=dbs.end()) return it->second;
    std::ifstream f(std::filesystem::path(MATERIALS_DATA_PATH) / file);
    if (!f) throw std::runtime_error(std::string("Materials: cannot open ") + file + " under " MATERIALS_DATA_PATH);
    json j; f >> j;
    return dbs.emplace(file, std::move(j)).first->second;
}

std::vector<std::string> EntryNames(const json& db)
{
    std::vector<std::string> names;
    for (const auto& [k,v] : db.items()) if (!k.empty() && k[0]!='_') names.push_back(k);
    return names;
}

const json& Entry(const char* file, const std::string& name)
{
    const json& db=Database(file);
    auto it=db.find(name);
    if (it==db.end() || name.empty() || name[0]=='_')
    {
        std::string known;
        for (const auto& n : EntryNames(db)) known += (known.empty() ? "" : ", ") + n;
        throw std::runtime_error("Materials: no entry '" + name + "' in " + file + "; known: " + known);
    }
    return *it;
}

int ZOf(const std::string& el)
{
    const size_t Z=thePeriodicTable().GetZ(el);
    if (Z==0) throw std::runtime_error("Materials: unknown element symbol '" + el + "'");
    return int(Z);
}

Bravais TypeOf(const std::string& s)
{
    static const std::map<std::string,Bravais> types = {
        {"CubicP",Bravais::CubicP}, {"CubicI",Bravais::CubicI}, {"CubicF",Bravais::CubicF},
        {"TetragonalP",Bravais::TetragonalP}, {"TetragonalI",Bravais::TetragonalI},
        {"OrthorhombicP",Bravais::OrthorhombicP}, {"OrthorhombicC",Bravais::OrthorhombicC},
        {"OrthorhombicI",Bravais::OrthorhombicI}, {"OrthorhombicF",Bravais::OrthorhombicF},
        {"HexagonalP",Bravais::HexagonalP}, {"RhombohedralR",Bravais::RhombohedralR},
        {"MonoclinicP",Bravais::MonoclinicP}, {"MonoclinicC",Bravais::MonoclinicC}, {"TriclinicP",Bravais::TriclinicP}};
    auto it=types.find(s);
    if (it==types.end()) throw std::runtime_error("Materials: unknown Bravais type '" + s + "'");
    return it->second;
}

Material FromEntry(const std::string& name, const json& e, double aOverride)
{
    const json& lat=e.at("lattice");
    LatticeParams p;
    p.a = aOverride>0.0 ? aOverride : lat.at("a").get<double>();
    if (lat.contains("b")) p.b=lat["b"].get<double>();
    if (lat.contains("c")) p.c=lat["c"].get<double>();
    if (lat.contains("alpha")) p.α=lat["alpha"].get<double>();
    if (lat.contains("beta"))  p.β=lat["beta"] .get<double>();
    if (lat.contains("gamma")) p.γ=lat["gamma"].get<double>();
    Matrix3D<int> T;
    if (lat.contains("T"))
    {
        const json& t=lat["T"];
        if (t.size()!=3 || t[0].size()!=3) throw std::runtime_error("Materials: '" + name + "': T must be a 3x3 integer matrix");
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) T(i+1,j+1)=t[i][j].get<int>();
    }
    Material m;
    m.name=name;
    m.cell=std::make_shared<UnitCell>(BravaisCell(TypeOf(lat.at("type").get<std::string>()), p, T));
    for (const auto& a : e.at("atoms"))
    {
        const json& f=a.at("frac");
        const int spin = a.contains("spin") ? a["spin"].get<int>() : 0;
        if (spin!=0 && spin!=1 && spin!=-1) throw std::runtime_error("Materials: '" + name + "': a site spin is +1, -1 or absent");
        m.cell->AddAtom(ZOf(a.at("el").get<std::string>()), rvec3_t(f[0].get<double>(), f[1].get<double>(), f[2].get<double>()),
                        /*spinFlip*/ spin==-1);
    }
    for (const auto& [el,val] : e.at("species").items()) m.species.emplace_back(el, val.get<int>());
    return m;
}
} // anonymous

int Material::Nelec() const
{
    int n=0;
    cell->ForEachSite([&](int Z, const rvec3_t&, bool)
    {
        for (const auto& [el,val] : species)
            if (ZOf(el)==Z) { n+=val; return; }
        throw std::runtime_error("Materials: '" + name + "': atom Z=" + std::to_string(Z) + " has no species entry");
    });
    return n;
}

Material Get(const std::string& name, double a) { return FromEntry(name, Entry("materials.json", name), a); }
std::vector<std::string> Names() { return EntryNames(Database("materials.json")); }

Material AtomInBox(const std::string& element, int valence, double a)
{
    Material m;
    m.name=element+"_box";
    m.cell=std::make_shared<UnitCell>(a);
    m.cell->AddAtom(ZOf(element), rvec3_t(0.5,0.5,0.5));
    m.species={{element, valence}};
    return m;
}

Material DimerInBox(const std::string& element, int valence, double a, double d, bool afm)
{
    Material m;
    m.name=element+"2_box";
    m.cell=std::make_shared<UnitCell>(a);
    m.cell->AddAtom(ZOf(element), rvec3_t(0.5-0.5*d/a,0.5,0.5), false);
    m.cell->AddAtom(ZOf(element), rvec3_t(0.5+0.5*d/a,0.5,0.5), afm);
    m.species={{element, valence}};
    return m;
}

Molecule GetMolecule(const std::string& name)
{
    const json& e=Entry("molecules.json", name);
    Molecule mol;
    for (const auto& a : e.at("atoms"))
    {
        const json& r=a.at("xyz");
        mol.Insert(new Atom(ZOf(a.at("el").get<std::string>()), 0.0, rvec3_t(r[0].get<double>(), r[1].get<double>(), r[2].get<double>())));
    }
    return mol;
}
std::vector<std::string> MoleculeNames() { return EntryNames(Database("molecules.json")); }

} // namespace qchem::Materials
