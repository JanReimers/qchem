// File: BasisSet/Lattice_3D/Imp/GTH_Potentials.C  GTH database reader implementation (JSON -> potentials).
module;
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Lattice_3D.GTH_Potentials;

#ifndef LATTICE_DATA_PATH
#error "LATTICE_DATA_PATH must be defined by CMake"
#endif

namespace BasisSet::Lattice_3D
{

//! The parsed database, loaded once on first use (construct-on-first-use; no static-init-order issue).
static const nlohmann::json& database()
{
    static const nlohmann::json db = []
    {
        std::ifstream f(std::filesystem::path(LATTICE_DATA_PATH) / "gth_potentials.json");
        if (!f) throw std::runtime_error("GTH: cannot open gth_potentials.json under " LATTICE_DATA_PATH);
        nlohmann::json j; f >> j; return j;
    }();
    return db;
}

GTH_PP GetGTH(const std::string& element, const std::string& functional, int q)
{
    const nlohmann::json& db = database();
    if (!db.contains(element))
        throw std::runtime_error("GTH: element '" + element + "' not in the database");
    const nlohmann::json& fns = db[element];
    if (!fns.contains(functional))
        throw std::runtime_error("GTH: " + element + " has no '" + functional + "' pseudopotential");
    const nlohmann::json& byq = fns[functional];

    std::string key = (q > 0) ? std::to_string(q) : byq.at("default").get<std::string>();
    if (!byq.contains(key))
        throw std::runtime_error("GTH: " + element + " " + functional + " has no q=" + key);
    const nlohmann::json& rec = byq[key];

    int zion = rec["z_ion"].get<int>();
    std::vector<double> c = rec["c"].get<std::vector<double>>();
    c.resize(4, 0.0);                       // pad to the analytic form's C1..C4
    HGH_LocalPotential local(zion, rec["r_loc"].get<double>(), c[0], c[1], c[2], c[3]);

    HGH_SeparablePotential nonlocal;
    for (const nlohmann::json& ch : rec["channels"])
        nonlocal.AddChannel(ch["l"].get<int>(), ch["r"].get<double>(),
                            ch["h"].get<std::vector<std::vector<double>>>());

    return GTH_PP{zion, local, nonlocal};
}

} // namespace
