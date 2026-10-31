// File: BasisSet/Gaussian/Point/Imp/Factory.C  Factory function for molecular basis sets (see Factory.C for the
// three orthogonal axes: BasisSetData x Engine x Angular).
module;
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <stdexcept>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Point.BasisFiles;
import qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Point.PG_Spherical;
import qchem.BasisSet.Gaussian.Point.PG_LibCint;

using json = nlohmann::json;

namespace qchem::BasisSet::Gaussian
{
    // The registries -- the single place each axis' names live.
    // THE ONE REGISTRATION of a C++ handle: enum -> file stem (the name everywhere else; the file is <stem>.bsd).
    static const std::map<BasisSetData, std::string> theStems =
    {
        {BasisSetData::DZVP , "dzvp" }, {BasisSetData::DZVP2, "dzvp2"}, {BasisSetData::TZVP, "tzvp"},
        {BasisSetData::ORB  , "orb"  }, {BasisSetData::ORB1 , "orb1" }, {BasisSetData::SIPP, "sipp"},
        {BasisSetData::SIPP_SR, "sipp_sr"},
        {BasisSetData::VALENCE_LOWQ,     "valence_lowq"    }, {BasisSetData::VALENCE_LOWQ_SR,  "valence_lowq_sr" },
        {BasisSetData::VALENCE_LOWQ_SR2, "valence_lowq_sr2"}, {BasisSetData::VALENCE_LOWQ_SPH, "valence_lowq_sph"},
        {BasisSetData::VALENCE_LOWQ_VA , "valence_lowq_va" }, {BasisSetData::VALENCE_LOWQ_VB , "valence_lowq_vb" },
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
            throw std::runtime_error(std::string("Gaussian::Factory: ")
                + (name.empty() ? std::string("missing required \"")+key+"\""
                                : std::string("unknown ")+key+" \""+name+"\"")
                + "; valid: " + valid);
        }
        return it->second;
    }

    // The declared {ELEMENT (upper-cased), q} of a file's `PP-PROVENANCE` preamble line; empty = nothing declared.
    // The preamble is everything before the first lone `!` (the Gaussian94 reader skips it, so the tag is invisible to it).
    static std::map<std::string,int> DeclaredPPProvenance(const std::string& path)
    {
        const auto upper=[](std::string s){ for (char& c : s) c=(char)std::toupper((unsigned char)c); return s; };
        std::ifstream f(path);
        if (!f) throw std::runtime_error("CheckPPProvenance: cannot open basis file " + path);
        std::map<std::string,int> declared;
        std::string line;
        while (std::getline(f,line))
        {
            std::istringstream is(line);
            std::string tok;
            if (!(is>>tok)) continue;
            if (tok=="!") break;                                   // end of the preamble
            if (tok!="PP-PROVENANCE") continue;
            while (is>>tok)
            {
                const size_t eq=tok.find('=');
                if (eq==std::string::npos || eq==0 || eq+1>=tok.size())
                    throw std::runtime_error("CheckPPProvenance: malformed PP-PROVENANCE entry '" + tok + "' in " + path + " (want Element=q)");
                declared[upper(tok.substr(0,eq))] = std::stoi(tok.substr(eq+1));
            }
        }
        return declared;
    }

    void CheckPPProvenance(const std::string& bsdFile, const std::vector<std::pair<std::string,int>>& species)
    {
        const std::string path=BasisFile(bsdFile);
        const std::map<std::string,int> declared=DeclaredPPProvenance(path);
        for (const auto& [element,q] : species)
        {
            std::string key=element; for (char& c : key) c=(char)std::toupper((unsigned char)c);
            const auto it=declared.find(key);
            if (it!=declared.end() && it->second!=q)
                throw std::runtime_error("basis/pseudopotential mismatch: " + bsdFile + " declares its " + element + " block for q"
                    + std::to_string(it->second) + " but the run's pseudopotential for " + element + " is q" + std::to_string(q)
                    + " -- regenerate the block with valgen, or choose the pseudopotential variant the file was made for");
        }
    }

    void CheckPPProvenance(BasisSetData data, const std::vector<std::pair<std::string,int>>& species)
    {
        CheckPPProvenance(BasisName(data)+".bsd", species);
    }

    std::string BasisName(BasisSetData data) {return theStems.at(data);}

    bool RegisteredBasis(const std::string& name, BasisSetData& data)
    {
        std::string stem=name;
        for (char& c : stem) c=(char)std::tolower((unsigned char)c);
        for (const auto& [d,s] : theStems) if (s==stem) { data=d; return true; }
        return false;
    }

    tBasisSet<double>* Factory(BasisSetData data, const Structure* cl, Engine engine, Angular angular,
                               const ShellTrim& trim)
    {
        return FactoryFromFile(BasisName(data)+".bsd", cl, engine, angular, trim);
    }

    tBasisSet<double>* FactoryFromFile(const std::string& bsdFile, const Structure* cl, Engine engine, Angular angular,
                                       const ShellTrim& trim)
    {
        Gaussian94Reader file(BasisFile(bsdFile));
        TrimmingReader reader(file, trim);   // pin 22: the vet-stage trim applied AS the file is read
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

    tBasisSet<double>* Factory(const json& js, const Structure* cl)
    {
        if (!js.contains("basis") || !js["basis"].is_string() || js["basis"].get<std::string>().empty())
            throw std::runtime_error("Gaussian::Factory: missing required \"basis\" (a file stem; the names: BasisNames())");
        const std::string file = BasisFileOfName(js["basis"].get<std::string>());   // required; any stem in the data directory
        Engine       engine  = Resolve(js, "engine",  "mnd",       theEngines);
        Angular      angular = Resolve(js, "angular", "cartesian", theAngulars);
        return FactoryFromFile(file, cl, engine, angular);
    }
}
