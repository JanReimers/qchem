// File: Structure/Imp/StructureData.C  The JSON reader behind qchem.StructureData (structure fields ONLY).
module;
#include <fstream>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <nlohmann/json.hpp>
module qchem.StructureData;
import qchem.PeriodicTable;   // thePeriodicTable().GetZ(symbol)
import qchem.Matrix3D;
import qchem.Types;

namespace qchem::StructureData
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
    std::ifstream f(std::filesystem::path(STRUCTURE_DATA_PATH) / file);
    if (!f) throw std::runtime_error(std::string("StructureData: cannot open ") + file + " under " STRUCTURE_DATA_PATH);
    json j; f >> j;
    return dbs.emplace(file, std::move(j)).first->second;
}

std::vector<std::string> EntryNames(const json& db)
{
    std::vector<std::string> names;
    for (const auto& [k,v] : db.items()) if (!k.empty() && k[0]!='_') names.push_back(k);
    return names;
}

bool Has(const json& db, const std::string& name) { return !name.empty() && name[0]!='_' && db.contains(name); }

std::string KnownNames()
{
    std::string known;
    for (const char* f : {"molecules.json","materials.json"})
        for (const auto& n : EntryNames(Database(f))) known += (known.empty() ? "" : ", ") + n;
    return known;
}

int ZOf(const std::string& el)
{
    const size_t Z=thePeriodicTable().GetZ(el);
    if (Z==0) throw std::runtime_error("StructureData: unknown element symbol '" + el + "'");
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
    if (it==types.end()) throw std::runtime_error("StructureData: unknown Bravais type '" + s + "'");
    return it->second;
}
} // anonymous

Kind KindOf(const std::string& name)
{
    const bool mol=Has(Database("molecules.json"), name), cell=Has(Database("materials.json"), name);
    if (mol && cell) throw std::runtime_error("StructureData: '" + name + "' is in BOTH molecules.json and materials.json");
    if (mol)  return Kind::Molecule;
    if (cell) return Kind::Cell;
    throw std::runtime_error("StructureData: no entry '" + name + "'; known: " + KnownNames());
}

Molecule GetMolecule(const std::string& name)
{
    if (KindOf(name)!=Kind::Molecule) throw std::runtime_error("StructureData: '" + name + "' is a periodic CELL, not a molecule (use GetCell)");
    const json& e=Database("molecules.json").at(name);
    Molecule mol;
    for (const auto& a : e.at("atoms"))
    {
        const json& r=a.at("xyz");
        mol.Insert(new Atom(ZOf(a.at("el").get<std::string>()), 0.0, rvec3_t(r[0].get<double>(), r[1].get<double>(), r[2].get<double>())));
    }
    return mol;
}
std::vector<std::string> MoleculeNames() { return EntryNames(Database("molecules.json")); }

UnitCell GetCell(const std::string& name, double aOverride)
{
    if (KindOf(name)!=Kind::Cell) throw std::runtime_error("StructureData: '" + name + "' is a MOLECULE, not a periodic cell (use GetMolecule)");
    const json& e=Database("materials.json").at(name);
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
        if (t.size()!=3 || t[0].size()!=3) throw std::runtime_error("StructureData: '" + name + "': T must be a 3x3 integer matrix");
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) T(i+1,j+1)=t[i][j].get<int>();
    }
    UnitCell cell=BravaisCell(TypeOf(lat.at("type").get<std::string>()), p, T);
    for (const auto& a : e.at("atoms"))
    {
        const json& f=a.at("frac");
        const int spin = a.contains("spin") ? a["spin"].get<int>() : 0;
        if (spin!=0 && spin!=1 && spin!=-1) throw std::runtime_error("StructureData: '" + name + "': a site spin is +1, -1 or absent");
        cell.AddAtom(ZOf(a.at("el").get<std::string>()), rvec3_t(f[0].get<double>(), f[1].get<double>(), f[2].get<double>()),
                     /*spinFlip*/ spin==-1);
    }
    return cell;
}
std::vector<std::string> CellNames() { return EntryNames(Database("materials.json")); }

UnitCell AtomInBox(const std::string& element, double a)
{
    UnitCell c(a);
    c.AddAtom(ZOf(element), rvec3_t(0.5,0.5,0.5));
    return c;
}

UnitCell DimerInBox(const std::string& element, double a, double d, bool afm)
{
    UnitCell c(a);
    c.AddAtom(ZOf(element), rvec3_t(0.5-0.5*d/a,0.5,0.5), false);
    c.AddAtom(ZOf(element), rvec3_t(0.5+0.5*d/a,0.5,0.5), afm);
    return c;
}

} // namespace qchem::StructureData
