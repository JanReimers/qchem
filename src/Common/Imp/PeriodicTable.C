// File: Common/Imp/PeriodicTable.C
module;
#include <string>
#include <cassert>
#include <vector> 
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <nlohmann/json.hpp>

module Common.PeriodicTable;


std::ostream& operator<<(std::ostream& os, const OrbitalRecordSaito& o)
{
    os << std::fixed;
    os << "    " << o.Symbol << " " << std::setw(14) << std::setprecision(11) << o.Energy_HF;
    os << std::fixed;
    for (auto r:o.r_moments)
        os << " " << std::setw(11) << std::setprecision(6) << r;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ElementRecordSaito& e)
{
    os << e.Symbol << "(" << e.Z << ") " << e.ValConfigString << " " << e.Term
    << " E_HF=" << std::setw(14) << std::setprecision(6) << e.Energy_HF
    << " E_DHF=" << std::setw(14) << std::setprecision(6) << e.Energy_DHF
    << " NUnpaired=" << e.NUnpaired << std::endl;
    // for (auto o:e.Orbitals)
    //     os << o << std::endl;
    return os;
}

OrbitalRecordSaito::OrbitalRecordSaito(nlohmann::json& j) : Symbol(j["name"]),Energy_HF(j["e"])
{
    for (auto r:j["rs"])
        r_moments.push_back(r);
}

std::set<size_t> fullShells({2,10,18,36,54,86,118});

// Relativistic orbital labels in physical (n,l,j) order, matching doc/DHF_GS_Energies_rel.json.
const std::vector<std::string> theDHFOrbitalLabels({
    "1s+",
    "2s+", "2p-", "2p+",
    "3s+", "3p-", "3p+", "3d-", "3d+",
    "4s+", "4p-", "4p+", "4d-", "4d+", "4f-", "4f+",
    "5s+", "5p-", "5p+", "5d-", "5d+", "5f-", "5f+",
    "6s+", "6p-", "6p+", "6d-", "6d+",
    "7s+",
});

#ifndef COMMON_DATA_PATH
#error "COMMON_DATA_PATH must be defined by CMake"
#endif

// Construct-on-first-use: this path is read from the ctors of static PeriodicTable* objects that may
// live in other translation units (e.g. QchemTester::itsPT).  A namespace-scope global would be a
// static-initialization-order fiasco -- it could still be zero-initialized when such a ctor runs,
// giving a null path string and a crash.  A function-local static is initialized on first call.
static const std::filesystem::path& common_data_dir()
{
    static const std::filesystem::path p = COMMON_DATA_PATH;
    return p;
}

ElementRecordSaito::ElementRecordSaito(nlohmann::json& j) : Z(j["Z"]), Symbol(j["symbol"]),  ValConfigString(j["valance"]), Term(j["term"]), MaxL(0), Energy_HF(j["HFEnergy"])
{
    for (auto o:j["Orbitals"])
        Orbitals.push_back(OrbitalRecordSaito(o));
    NUnpaired=std::atoi(&Term[0])-1;
    assert(NUnpaired<=8);
    size_t spos=ValConfigString.find('s'),ppos=ValConfigString.find('p'),dpos=ValConfigString.find('d'),fpos=ValConfigString.find('f');
    bool bs=spos!=std::string::npos,bp=ppos!=std::string::npos,bd=dpos!=std::string::npos,bf=fpos!=std::string::npos;
    std::string nss="0",nps="0",nds="0",nfs="0";

    if (bs)
        nss=ValConfigString.substr(spos+1,1); //number of s electron is always one digit.
    if (bp)
        nps=ValConfigString.substr(ppos+1,1); //number of p electron is always one digit.
    if (bd)
    {
        if (bp)
            nds=ValConfigString.substr(dpos+1,ppos-dpos-2); //d is always just before p.
        else
            nds=ValConfigString.substr(dpos+1,ValConfigString.length()-dpos-1);
    }
    if (bf)
    {
        if (bd)
            nfs=ValConfigString.substr(fpos+1,dpos-fpos-2); //f is usually just before d.
        else if (bp)
            nfs=ValConfigString.substr(fpos+1,ppos-fpos-2); //No d, so f is just before p.
        else
            nfs=ValConfigString.substr(fpos+1,ValConfigString.length()-fpos-1);
    }
    ValConfig[0]=stoi(nss);
    ValConfig[1]=stoi(nps);
    ValConfig[2]=stoi(nds);
    ValConfig[3]=stoi(nfs);

    if (Z>4) MaxL=1;
    if (Z>20) MaxL=2;
    if (Z>57) MaxL=3;

    for (size_t l:iv_t(0,4))
    {
        assert(ValConfig[l]<=2*(2*l+1));
        if (fullShells.find(Z)!=fullShells.end()) ValConfig[l]=0; //Clear out full shells.
    }
}


PeriodicTableSaito::PeriodicTableSaito()
{
    // Read in Saito HF data.
    {
        std::ifstream file(common_data_dir() / "saito.json");
        assert(file);
        nlohmann::json jsondata;
        file >> jsondata;
        for (auto e:jsondata)
            elements.push_back(ElementRecordSaito(e));
    }
    // Read in NIST DFT data
    {
        std::ifstream file(common_data_dir() / "nistLDA.json");
        assert(file);
        nlohmann::json jsondata;
        file >> jsondata;
        for (auto& jse:jsondata)
        {
            size_t Z=jse["Z"].template get<size_t>();
            ElementRecordSaito& e=elements[Z-1];
            e.Energy_DFT=jse["Etot"].template get<double>();
            for (auto& o:e.Orbitals)
            {
                o.Energy_DFT=jse[o.Symbol].template get<double>();
            }
        }
    }
    // Read in relativistic DHF data (total energy + spin-orbit split orbital eigenvalues).
    {
        std::ifstream file(common_data_dir() / "DHF_GS_Energies_rel.json");
        assert(file);
        nlohmann::json jsondata;
        file >> jsondata;
        for (auto& jse:jsondata)
        {
            size_t Z=jse["Z"].template get<size_t>();
            ElementRecordSaito& e=elements[Z-1];
            e.Energy_DHF=jse["TE_rel"].template get<double>();
            for (const std::string& label:theDHFOrbitalLabels)
            {
                const auto& val=jse[label];
                if (val.is_null()) continue;
                e.DHFOrbitals.push_back({label,val.template get<double>()});
            }
        }
    }

}

