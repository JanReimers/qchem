// File: BasisSet/Molecule/Imp/Factory.C  Factory function for molecular basis sets.
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
import qchem.BasisSet.Molecule.PG_Cart_LibCint;

using json = nlohmann::json;

namespace BasisSet::Molecule
{
    // The orbital-basis catalogue: enum -> data file, and name -> enum (for the json/config path).  The
    // single place a basis is registered.
    static const std::map<Basis, std::string> theFiles =
    {
        {Basis::DZVP , "dzvp.bsd" },
        {Basis::DZVP2, "dzvp2.bsd"},
        {Basis::TZVP , "tzvp.bsd" },
        {Basis::ORB  , "orb.bsd"  },
        {Basis::ORB1 , "orb1.bsd" },
    };
    static const std::map<std::string, Basis> theNames =
    {
        {"dzvp" , Basis::DZVP }, {"dzvp2", Basis::DZVP2}, {"tzvp", Basis::TZVP},
        {"orb"  , Basis::ORB  }, {"orb1" , Basis::ORB1 },
    };

    BasisSet<double>* Factory(Basis basis, const Cluster* cl, bool spherical)
    {
        Gaussian94Reader reader(BasisFile(theFiles.at(basis)));
        if (spherical) return new PG_Spherical::BasisSet(&reader,cl);
        return new PG_Cart::BasisSet(&reader,cl);
    }

    BasisSet<double>* Factory(const nlohmann::json& js,const Cluster* cl)
    {
        const std::string name = js["basis"].template get<std::string>();
        auto it = theNames.find(name);
        if (it==theNames.end())
        {
            std::string valid;
            for (const auto& [n,b] : theNames) valid += (valid.empty()?"":", ") + n;
            throw std::runtime_error("Molecule::Factory: unknown basis \"" + name + "\"; valid: " + valid);
        }
        // Optional js["engine"]: "mnd" (default, McMurchie-Davidson) or "libcint" (Cartesian only, HF).
        if (js.value("engine", std::string("mnd")) == "libcint")
        {
            Gaussian94Reader reader(BasisFile(theFiles.at(it->second)));
            return new PG_Cart_LibCint::BasisSet(&reader, cl);
        }
        return Factory(it->second, cl, js.value("spherical", false));
    }
}
