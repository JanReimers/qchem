// File: PeriodicTable.C  Test the periodic table class, or database.
#include "gtest/gtest.h"
#include <vector>
#include <fstream>
#include <iomanip>

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