// File: Structure/Imp/StructureData.C  The JSON reader behind qchem.StructureData (structure fields ONLY).
module;
#include <algorithm>
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
    throw std::runtime_error("StructureData: no entry '" + name + "' among the shipped structures; known: " + KnownNames()
                             + " -- a user material must be inlined in the deck (its entry as the structure), never named");
}

Molecule GetMolecule(const std::string& name)
{
    if (KindOf(name)!=Kind::Molecule) throw std::runtime_error("StructureData: '" + name + "' is a periodic CELL, not a molecule (use GetCell)");
    const json& e=Database("molecules.json").at(name);
    Molecule mol;
    for (const auto& a : e.at("atoms"))
    {
        const json& r=a.at("xyz");
        mol.Insert(Site{ZOf(a.at("el").get<std::string>()), rvec3_t(r[0].get<double>(), r[1].get<double>(), r[2].get<double>())});
    }
    return mol;
}
std::vector<std::string> MoleculeNames() { return EntryNames(Database("molecules.json")); }

namespace
{
// ONE parser for a materials.json-schema entry (CellOf below), the data file's and an inline one's (a deck's object `structure`).
// Strict keys: a typo ("alfa", "fracs") would otherwise build a DIFFERENT cell silently.
void OnlyKeys(const json& o, std::initializer_list<const char*> keys, const std::string& where, const std::string& name)
{
    if (!o.is_object()) throw std::runtime_error("StructureData: '" + name + "': " + where + " must be an object");
    for (const auto& [k,v] : o.items())
        if (std::none_of(keys.begin(), keys.end(), [&](const char* x){ return k==x; }))
        {
            std::string known; for (const char* x : keys) known += (known.empty() ? "" : ", ") + std::string(x);
            throw std::runtime_error("StructureData: '" + name + "': unknown key '" + k + "' in " + where + " (known: " + known + ")");
        }
}

// The atom basis of an entry onto \a cell (fractional coordinates of THAT cell, an optional collinear spin flip).
UnitCell AddAtoms(UnitCell cell, const json& e, const std::string& name)
{
    for (const auto& a : e.at("atoms"))
    {
        OnlyKeys(a, {"el","frac","spin"}, "an atom", name);
        const json& f=a.at("frac");
        if (!f.is_array() || f.size()!=3) throw std::runtime_error("StructureData: '" + name + "': an atom's 'frac' must be [f1,f2,f3]");
        const int spin = a.contains("spin") ? a["spin"].get<int>() : 0;
        if (spin!=0 && spin!=1 && spin!=-1) throw std::runtime_error("StructureData: '" + name + "': a site spin is +1, -1 or absent");
        cell.AddAtom(ZOf(a.at("el").get<std::string>()), rvec3_t(f[0].get<double>(), f[1].get<double>(), f[2].get<double>()),
                     /*spinFlip*/ spin==-1);
    }
    return cell;
}

UnitCell CellOf(const json& e, const std::string& name, double aOverride)
{
    auto only=[&](const json& o, std::initializer_list<const char*> keys, const std::string& where){ OnlyKeys(o, keys, where, name); };
    if (!e.contains("lattice") || !e.contains("atoms") || !e.at("atoms").is_array() || e.at("atoms").empty())
        throw std::runtime_error("StructureData: '" + name + "': an entry needs 'lattice' and a non-empty 'atoms' list");
    const json& lat=e.at("lattice");
    // EXPLICIT VECTORS (2026-10-31, the GUI's inline cells): `vectors_bohr` [[a],[b],[c]], the three lattice vectors in bohr --
    // exclusive with the Bravais description (type + parameters + T), which is the data file's own shorthand for the same thing.
    if (lat.contains("vectors_bohr"))
    {
        only(lat, {"vectors_bohr"}, "lattice (vectors_bohr is exclusive with type / a..gamma / T)");
        const json& v=lat.at("vectors_bohr");
        if (!v.is_array() || v.size()!=3) throw std::runtime_error("StructureData: '" + name + "': vectors_bohr must be [[a],[b],[c]], three 3-vectors");
        Matrix3D<double> A;
        for (int i=0;i<3;i++)
        {
            if (!v[i].is_array() || v[i].size()!=3) throw std::runtime_error("StructureData: '" + name + "': vectors_bohr must be [[a],[b],[c]], three 3-vectors");
            for (int x=0;x<3;x++) A(x+1,i+1)=v[i][x].get<double>();   // the cell matrix's COLUMNS are the lattice vectors
        }
        if (aOverride>0.0) throw std::runtime_error("StructureData: '" + name + "': a lattice-constant override needs a Bravais lattice, not vectors_bohr");
        return AddAtoms(UnitCell(A), e, name);
    }
    only(lat, {"type","a","b","c","alpha","beta","gamma","T"}, "lattice");
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
    return AddAtoms(BravaisCell(TypeOf(lat.at("type").get<std::string>()), p, T), e, name);
}
} // anonymous

UnitCell GetCell(const std::string& name, double aOverride)
{
    if (KindOf(name)!=Kind::Cell) throw std::runtime_error("StructureData: '" + name + "' is a MOLECULE, not a periodic cell (use GetMolecule)");
    return CellOf(Database("materials.json").at(name), name, aOverride);
}

UnitCell CellFromEntry(const std::string& entryJson, const std::string& who)
{
    json e;
    try { e=json::parse(entryJson); }
    catch (const json::exception& x) { throw std::runtime_error("StructureData: '" + who + "': not JSON: " + x.what()); }
    return CellOf(e, who, 0.0);
}
std::vector<std::string> CellNames() { return EntryNames(Database("materials.json")); }

// ---------------------------------------------------------------------------------------------- THE LAYERS (OpenWork §2, 2026-11-09)
namespace
{
//! A library file (materials.json schema), re-read when it changed on disk; a missing file is an empty layer.
const json& Library(const std::string& path)
{
    static std::map<std::string, std::pair<std::filesystem::file_time_type, json>> cache;
    static const json none=json::object();
    std::error_code ec;
    const auto t=std::filesystem::last_write_time(path, ec);
    if (ec) { cache.erase(path); return none; }
    if (auto it=cache.find(path); it!=cache.end() && it->second.first==t) return it->second.second;
    std::ifstream f(path);
    json j;
    try { f >> j; }
    catch (const json::exception& x) { throw std::runtime_error("StructureData: materials library " + path + " is not JSON: " + x.what()); }
    if (!j.is_object()) throw std::runtime_error("StructureData: materials library " + path + " must be an object of named entries");
    return (cache[path]={t, std::move(j)}).second;
}

std::string LayerName(Layer l) { return l==Layer::Shipped ? "the shipped materials" : l==Layer::User ? "your materials library" : "a workspace library"; }

void RequireLegalName(const std::string& name, const std::string& where)
{
    if (name.empty() || name[0]=='_' || name.find('/')!=std::string::npos)
        throw std::runtime_error("StructureData: " + where + ": name '" + name + "' must be non-empty, not start with '_' and not contain '/'");
}
} // anonymous

std::vector<CellEntry> CellEntries(const std::vector<LibraryFile>& layers)
{
    std::vector<CellEntry> out;
    const auto molecules=EntryNames(Database("molecules.json"));
    const std::string shipped=(std::filesystem::path(STRUCTURE_DATA_PATH) / "materials.json").string();
    for (const auto& n : EntryNames(Database("materials.json")))
        out.push_back({n, Layer::Shipped, shipped, Database("materials.json").at(n).dump()});
    for (const auto& [layer, path] : layers)
    {
        if (layer==Layer::Shipped) throw std::runtime_error("StructureData: " + path + ": a library layered on top is a user or workspace one, not shipped");
        const json& db=Library(path);
        for (const auto& n : EntryNames(db))
        {
            if (std::find(molecules.begin(), molecules.end(), n)!=molecules.end())
                throw std::runtime_error("StructureData: " + path + ": '" + n + "' is a shipped molecule name -- rename it in " + path);
            for (auto& e : out) if (e.name==n) e.shadowed=true;      // this later layer's entry is the one the name finds
            out.push_back({n, layer, path, db.at(n).dump()});
        }
    }
    return out;
}

CellEntry GetCellEntry(const std::string& name, const std::vector<LibraryFile>& layers)
{
    const auto all=CellEntries(layers);
    for (const auto& e : all) if (e.name==name && !e.shadowed) return e;
    std::string known; for (const auto& e : all) if (!e.shadowed) known += (known.empty() ? "" : ", ") + e.name;
    throw std::runtime_error("StructureData: no material '" + name + "'; known: " + known
                             + (layers.empty() ? " -- a user material must be inlined in the deck (its entry as the structure), never named" : ""));
}

void SaveMaterial(const LibraryFile& into, const std::string& name, const std::string& entryJson, bool replace)
{
    if (into.layer==Layer::Shipped)
        throw std::runtime_error("StructureData: the shipped materials are CURATED -- promoting a material into them is a deliberate commit, not a save");
    RequireLegalName(name, into.path);
    json e;
    try { e=json::parse(entryJson); }
    catch (const json::exception& x) { throw std::runtime_error("StructureData: '" + name + "' is not JSON: " + x.what()); }
    if (!e.is_object()) throw std::runtime_error("StructureData: '" + name + "' must be an object (a materials.json entry)");
    e.erase("name");                                  // the KEY names it
    (void)CellOf(e, name, 0.0);                       // the parser refuses before the file is touched
    if (Has(Database("molecules.json"), name)) throw std::runtime_error("StructureData: '" + name + "' is a shipped molecule name -- choose another");
    json db=Library(into.path);
    if (db.contains(name) && !replace)
        throw std::runtime_error("StructureData: '" + name + "' is already in " + LayerName(into.layer) + " -- replace it, or choose another name");
    db[name]=e;                                       // ordered: a replaced entry keeps its place, a new one goes last
    const std::filesystem::path p(into.path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
    const std::filesystem::path tmp=p.string()+".tmp";
    {
        std::ofstream f(tmp);
        if (!f) throw std::runtime_error("StructureData: cannot write " + tmp.string());
        f << db.dump(2) << "\n";
        if (!f) throw std::runtime_error("StructureData: writing " + tmp.string() + " failed");
    }
    std::filesystem::rename(tmp, p);                  // whole or not at all
}

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
