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

struct ElementRecord
    {
        size_t      Z;
        std::string Symbol;
        size_t      NUnpaired;
        size_t      MaxL;
        size_t      ValConfig[4]; //spdf
        double      EnergyHF;     //Saito, Shiro L. Hartree–Fock–Roothaan energies and expectation values for the neutral atoms He to Uuo: The B-spline expansion method, Atomic Data and Nuclear Data Tables, 95,6, 836--870
        double      EnergyDFT;    //NIST https://math.nist.gov/DFTdata/atomdata/tables/ptable.html
    };

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


TEST_F(PeriodicTableTests,ReadSaito)
{
    std::ifstream file("../../../src/Common/Data/saito.json");
    assert(file);
    nlohmann::json jsondata;
    file >> jsondata;
    std::vector<ElementRecordSaito> pt;
    for (auto e:jsondata)
        pt.push_back(ElementRecordSaito(e));
    for (auto e:pt)
        std::cout << e;
    
}

TEST_F(PeriodicTableTests,PeriodicTableSaito)
{
    PeriodicTableSaito pts;
    PeriodicTable pt;
    for (size_t Z:iv_t(1,92+1))
    {
        std::cout << Z << " " << pts.GetSymbol(Z) << std::endl;
        EXPECT_EQ(pt.GetSymbol(Z),pts.GetSymbol(Z));
        EXPECT_EQ(pt.GetEnergyHF(Z),pts.GetEnergyHF(Z));
        if (Z!=58) //The ground state configuration for Ce is controversial.
            EXPECT_EQ(pt.GetNumUnpairedElectrons(Z),pts.GetNumUnpairedElectrons(Z));
        EXPECT_EQ(pt.GetMaxL(Z),pts.GetMaxL(Z));
        for (size_t l:iv_t(0,3+1))
        {
            EXPECT_EQ(pt.GetValanceConfiguration(Z)[l],pts.GetValanceConfiguration(Z)[l]);
        }
    }
}