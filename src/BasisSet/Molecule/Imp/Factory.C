// File: BasisSet/Molecule/Imp/Factory.C  Factory function for molecular basis sets (see Factory.C for the
// three orthogonal axes: BasisSetData x Engine x Angular).
module;
#include <map>
#include <string>
#include <stdexcept>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Molecule.Factory;
import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.BasisFiles;
import qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet.Molecule.PG_Spherical;
import qchem.BasisSet.Molecule.PG_LibCint;

using json = nlohmann::json;

namespace BasisSet::Molecule
{
    // The registries -- the single place each axis' names live.
    static const std::map<BasisSetData, std::string> theFiles =      // axis 1: enum -> data file
    {
        {BasisSetData::DZVP , "dzvp.bsd" },
        {BasisSetData::DZVP2, "dzvp2.bsd"},
        {BasisSetData::TZVP , "tzvp.bsd" },
        {BasisSetData::ORB  , "orb.bsd"  },
        {BasisSetData::ORB1 , "orb1.bsd" },
    };
    static const std::map<std::string, BasisSetData> theBasisNames = // axis 1: json name -> enum
    {
        {"dzvp" , BasisSetData::DZVP }, {"dzvp2", BasisSetData::DZVP2}, {"tzvp", BasisSetData::TZVP},
        {"orb"  , BasisSetData::ORB  }, {"orb1" , BasisSetData::ORB1 },
    };
    static const std::map<std::string, Engine>  theEngines  =        // axis 2: json name -> enum
    { {"mnd", Engine::MnD}, {"libcint", Engine::LibCint} };
    static const std::map<std::string, Angular> theAngulars =        // axis 3: json name -> enum
    { {"cartesian", Angular::Cartesian}, {"spherical", Angular::Spherical} };

    // Resolve a json string key against a name->enum registry, throwing a clear error listing the valid
    // values.  `dflt` is used when the key is absent ("" = the key is required: an empty/absent value fails).
    template <class E>
    static E Resolve(const json& js, const char* key, const std::string& dflt,
                     const std::map<std::string,E>& names)
    {
        const std::string name = js.value(key, dflt);
        auto it = names.find(name);
        if (it==names.end())
        {
            std::string valid; for (const auto& [n,e]:names) valid += (valid.empty()?"":", ") + n;
            throw std::runtime_error(std::string("Molecule::Factory: ")
                + (name.empty() ? std::string("missing required \"")+key+"\""
                                : std::string("unknown ")+key+" \""+name+"\"")
                + "; valid: " + valid);
        }
        return it->second;
    }

    BasisSet<double>* Factory(BasisSetData data, const Structure* cl, Engine engine, Angular angular)
    {
        Gaussian94Reader reader(BasisFile(theFiles.at(data)));
        const bool spherical = (angular==Angular::Spherical);
        switch (engine)
        {
        case Engine::LibCint: return new PG_LibCint::BasisSet(&reader, cl, spherical);
        case Engine::MnD:     break;
        }
        // MnD: one IBS tree per angular kind.
        if (spherical) return new PG_Spherical::BasisSet(&reader, cl);
        return new PG_Cart::BasisSet(&reader, cl);
    }

    BasisSet<double>* Factory(const json& js, const Structure* cl)
    {
        BasisSetData data    = Resolve(js, "basis",   "",          theBasisNames);  // required
        Engine       engine  = Resolve(js, "engine",  "mnd",       theEngines);
        Angular      angular = Resolve(js, "angular", "cartesian", theAngulars);
        return Factory(data, cl, engine, angular);
    }
}
