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

// The PSEUDOPOTENTIAL half of an entry (valence counts).  INTERIM: still read from materials.json's `species`
// field; D-STRUCTDATA step 2 moves it to its own pseudopotentials.json.  The structure half is qcStructure's.
Material FromEntry(const std::string& name, const json& e, double aOverride)
{
    Material m;
    m.name=name;
    m.cell=std::make_shared<UnitCell>(StructureData::GetCell(name, aOverride));
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
