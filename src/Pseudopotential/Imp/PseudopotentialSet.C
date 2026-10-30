// File: Pseudopotential/Imp/PseudopotentialSet.C  The reader behind qchem.Pseudopotential.Set.
module;
#include <fstream>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
module qchem.Pseudopotential.Set;

#ifndef PSEUDO_DATA_PATH
#error "PSEUDO_DATA_PATH must be defined by CMake"
#endif

namespace qchem::Pseudopotential
{

namespace
{
using json = nlohmann::ordered_json;   // ORDERED: file order is the pick-list order

const json& Database()
{
    static const json db = []
    {
        std::ifstream f(std::filesystem::path(PSEUDO_DATA_PATH) / "pseudopotentials.json");
        if (!f) throw std::runtime_error("PseudopotentialSet: cannot open pseudopotentials.json under " PSEUDO_DATA_PATH);
        json j; f >> j; return j;
    }();
    return db;
}
} // anonymous

struct SetReader
{
    static PseudopotentialSet Read(const std::string& name, const json& e)
    {
        std::vector<std::pair<std::string,int>> valence;
        for (const auto& [el,q] : e.at("valence").items()) valence.emplace_back(el, q.get<int>());
        return PseudopotentialSet(name, e.at("functional").get<std::string>(), std::move(valence));
    }
};

std::vector<std::string> PseudopotentialSet::Names()
{
    std::vector<std::string> names;
    for (const auto& [k,v] : Database().items()) if (!k.empty() && k[0]!='_') names.push_back(k);
    return names;
}

PseudopotentialSet PseudopotentialSet::Get(const std::string& name)
{
    const json& db=Database();
    auto it=db.find(name);
    if (it==db.end() || name.empty() || name[0]=='_')
    {
        std::string known;
        for (const auto& n : Names()) known += (known.empty() ? "" : ", ") + n;
        throw std::runtime_error("PseudopotentialSet: no set '" + name + "' in pseudopotentials.json; known: " + known);
    }
    return SetReader::Read(name, *it);
}

bool PseudopotentialSet::Covers(const std::string& element) const
{
    for (const auto& [el,q] : itsValence) if (el==element) return true;
    return false;
}

int PseudopotentialSet::Valence(const std::string& element) const
{
    for (const auto& [el,q] : itsValence) if (el==element) return q;
    throw std::runtime_error("PseudopotentialSet '" + itsName + "' has no entry for element '" + element + "'");
}

std::vector<std::pair<std::string,int>> PseudopotentialSet::Species(const std::vector<std::string>& elements) const
{
    std::vector<std::pair<std::string,int>> s;
    for (const auto& el : elements) s.emplace_back(el, Valence(el));
    return s;
}

} // namespace qchem::Pseudopotential
