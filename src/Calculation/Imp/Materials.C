// File: Calculation/Imp/Materials.C  The JSON readers behind qchem.Materials.
module;
#include <algorithm>
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
import qchem.Pseudopotential.Set;   // the pseudopotential vocabulary a material names
import qchem.StructureData;   // GetCell/AtomInBox/DimerInBox -- the structure half of an entry

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
    std::ifstream f(std::filesystem::path(STRUCTURE_DATA_PATH) / file);
    if (!f) throw std::runtime_error(std::string("Materials: cannot open ") + file + " under " STRUCTURE_DATA_PATH);
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

// The PSEUDOPOTENTIAL half of an entry: the NAME of a set (`pp`), resolved against the cell's elements -- PP data is a property of
// (element, pseudopotential), not of the material (D-STRUCTDATA step 2).  The structure half is qcStructure's.
Material FromCell(const std::string& name, UnitCell cell, const std::string& pp)
{
    Material m;
    m.name=name;
    m.cell=std::make_shared<UnitCell>(std::move(cell));
    m.pp=pp;
    std::vector<std::string> elements;   // in order of first appearance in the cell
    m.cell->ForEachSite([&](int Z, const rvec3_t&, bool)
    {
        const std::string el=thePeriodicTable().GetSymbol(Z);
        if (std::find(elements.begin(), elements.end(), el)==elements.end()) elements.push_back(el);
    });
    m.species=Pseudopotential::PseudopotentialSet::Get(m.pp).Species(elements);
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

namespace
{
Material FromEntry(const std::string& name, const json& e, double aOverride)
{
    return FromCell(name, StructureData::GetCell(name, aOverride), e.at("pp").get<std::string>());
}
} // anonymous

Material Get(const std::string& name, double a) { return FromEntry(name, Entry("materials.json", name), a); }

Material Inline(const std::string& entryJson, const std::string& pp)
{
    json e;
    try { e=json::parse(entryJson); }
    catch (const json::exception& x) { throw std::runtime_error(std::string("Materials: an inline structure is not JSON: ") + x.what()); }
    if (!e.is_object()) throw std::runtime_error("Materials: an inline structure must be an object (a materials.json entry)");
    for (const auto& [k,v] : e.items())
        if (k!="name" && k!="comment" && k!="lattice" && k!="atoms" && k!="pp" && k!="provenance")
            throw std::runtime_error("Materials: unknown key '" + k + "' in an inline structure (known: name, comment, lattice, atoms, pp, provenance)");
    const std::string name = e.contains("name") ? e.at("name").get<std::string>() : std::string("inline");
    if (name.empty() || name[0]=='_' || name.find('/')!=std::string::npos)
        throw std::runtime_error("Materials: inline structure name '" + name + "' must be non-empty, not start with '_' and not contain '/' (it is a file stem)");
    const std::string set = e.contains("pp") ? e.at("pp").get<std::string>() : pp;
    if (set.empty()) throw std::runtime_error("Materials: inline structure '" + name + "' names no pseudopotential set (its 'pp', or the caller's)");
    return FromCell(name, StructureData::CellFromEntry(entryJson, name), set);
}
std::vector<std::string> Names() { return StructureData::CellNames(); }

Material AtomInBox(const std::string& element, int valence, double a)
{
    Material m;
    m.name=element+"_box";
    m.cell=std::make_shared<UnitCell>(StructureData::AtomInBox(element, a));
    m.species={{element, valence}};
    return m;
}

Material DimerInBox(const std::string& element, int valence, double a, double d, bool afm)
{
    Material m;
    m.name=element+"2_box";
    m.cell=std::make_shared<UnitCell>(StructureData::DimerInBox(element, a, d, afm));
    m.species={{element, valence}};
    return m;
}

} // namespace qchem::Materials
