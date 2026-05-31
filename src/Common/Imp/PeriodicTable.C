// File: Common/Imp/PeriodicTable.C
module;
#include <string>
#include <cassert>
#include <vector> 
#include <fstream>
#include <set>
#include <nlohmann/json.hpp>

module Common.PeriodicTable;


std::ostream& operator<<(std::ostream& os, const OrbitalRecordSaito& o)
{
    os << std::fixed;
    os << "    " << o.Symbol << " " << std::setw(14) << std::setprecision(11) << o.Energy;
    os << std::fixed;
    for (auto r:o.r_moments)
        os << " " << std::setw(11) << std::setprecision(6) << r;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ElementRecordSaito& e)
{
    os << e.Symbol << "(" << e.Z << ") " << e.ValConfigString << " " << e.Term 
    << " E_HF=" << std::setw(14) << std::setprecision(6) << e.EnergyHF 
    << " NUnpaired=" << e.NUnpaired << std::endl;
    // for (auto o:e.Orbitals)
    //     os << o << std::endl;
    return os;
}

OrbitalRecordSaito::OrbitalRecordSaito(nlohmann::json& j) : Symbol(j["name"]),Energy(j["e"])
{
    for (auto r:j["rs"])
        r_moments.push_back(r);
}

std::set<size_t> fullShells({2,10,18,36,54,86,118});

ElementRecordSaito::ElementRecordSaito(nlohmann::json& j) : Z(j["Z"]), Symbol(j["symbol"]),  ValConfigString(j["valance"]), Term(j["term"]), MaxL(0), EnergyHF(j["HFEnergy"])
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
    std::ifstream file("../../../doc/saito.json");
    assert(file);
    nlohmann::json jsondata;
    file >> jsondata;
    for (auto e:jsondata)
        elements.push_back(ElementRecordSaito(e));
}