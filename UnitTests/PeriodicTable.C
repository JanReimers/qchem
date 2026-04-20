// File: PeriodicTable.C  Test the periodic table class, or database.
#include "gtest/gtest.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

import Common.PeriodicTable;
import qchem.Types;

class PeriodicTableTests : public ::testing::Test
{};

TEST_F(PeriodicTableTests,DumpToCSV)
{
    PeriodicTable pt;
    std::vector<ElementRecord> recs;
    for (size_t Z:iv_t(1,92+1))
    {
        size_t NUnp=static_cast<size_t>(pt.GetNumUnpairedElectrons(Z));
        size_t MaxL=static_cast<size_t>(pt.GetMaxL(Z));
        auto vc=pt.GetValanceConfiguration(Z);
        ElementRecord rec={Z,pt.GetSymbol(Z),NUnp,MaxL
            ,{static_cast<size_t>(vc[0])
            ,static_cast<size_t>(vc[1])
            ,static_cast<size_t>(vc[2])
            ,static_cast<size_t>(vc[3])}
            ,pt.GetEnergyHF(Z),pt.GetEnergyDFT(Z)};
        recs.push_back(rec);
    }

    std::ofstream of("PeriodicTable.csv"); 
    of << "Z,Symbol,NUnpaired,MaxL,s,p,d,f,E_HF,E_DFT" << std::endl;
    for (auto rec:recs)
    {
        of << rec.Z << "," << rec.Symbol << "," << rec.NUnpaired << "," << rec.MaxL << ",";
        for (auto v:rec.ValConfig) of << v << ",";
        of << std::setprecision(13) << rec.EnergyHF << "," << rec.EnergyDFT;
        of << std::endl;
    }
}

struct OrbitalRecord
{
    OrbitalRecord(nlohmann::json& j) : Symbol(j["name"]),Energy(j["e"])
    {
        for (auto r:j["rs"])
            r_moments.push_back(r);
    }
    std::string Symbol;
    double      Energy;
    std::vector<double> r_moments; //<r^2>, <r^1>, <r^-1>, <r^-2>, <r^-3>, 
};

struct ElementRecord1
{
    ElementRecord1(nlohmann::json& j) : Z(j["Z"]), Symbol(j["symbol"]),  ValConfigString(j["valance"]), Term(j["term"]), EnergyHF(j["HFEnergy"])
    {
        for (auto o:j["Orbitals"])
            Orbitals.push_back(OrbitalRecord(o));
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
        for (size_t l:iv_t(0,4))
            assert(ValConfig[l]<=2*(2*l+1));
    }

    size_t      Z;
    std::string Symbol;
    std::string ValConfigString;
    std::string Term;
    size_t      NUnpaired;
    size_t      ValConfig[4]; //spdf
    double      EnergyHF;     //Saito, Shiro L. Hartree–Fock–Roothaan energies and expectation values for the neutral atoms He to Uuo: The B-spline expansion method, Atomic Data and Nuclear Data Tables, 95,6, 836--870
    std::vector<OrbitalRecord> Orbitals;
};
std::ostream& operator<<(std::ostream& os, const OrbitalRecord& o)
{
    os << std::fixed;
    os << "    " << o.Symbol << " " << std::setw(14) << std::setprecision(11) << o.Energy;
    os << std::fixed;
    for (auto r:o.r_moments)
        os << " " << std::setw(11) << std::setprecision(6) << r;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ElementRecord1& e)
{
    os << e.Symbol << "(" << e.Z << ") " << e.ValConfigString << " " << e.Term 
    << " E_HF=" << std::setw(14) << std::setprecision(6) << e.EnergyHF 
    << " NUnpaired=" << e.NUnpaired << std::endl;
    // for (auto o:e.Orbitals)
    //     os << o << std::endl;
    return os;
}

TEST_F(PeriodicTableTests,ReadSaito)
{
    std::ifstream file("../../../doc/saito.json");
    assert(file);
    nlohmann::json jsondata;
    file >> jsondata;
    std::vector<ElementRecord1> pt;
    for (auto e:jsondata)
        pt.push_back(ElementRecord1(e));
    for (auto e:pt)
        std::cout << e;
    
}