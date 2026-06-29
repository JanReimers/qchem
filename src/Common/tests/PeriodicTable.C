// File: PeriodicTable.C  Test the periodic table class, or database.
#include "gtest/gtest.h"
#include <vector>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>

import qchem.PeriodicTable;
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
    PeriodicTableSaito pt;
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


TEST_F(PeriodicTableTests,Electronegativity)
{
    const PeriodicTableSaito& pt = thePeriodicTable();
    // Spot-check the Z-indexing alignment across the table (a shift would misplace these).
    EXPECT_NEAR(pt.GetElectronegativity(1),  2.20, 1e-9);   // H
    EXPECT_NEAR(pt.GetElectronegativity(9),  3.98, 1e-9);   // F  (most electronegative)
    EXPECT_NEAR(pt.GetElectronegativity(11), 0.93, 1e-9);   // Na
    EXPECT_NEAR(pt.GetElectronegativity(17), 3.16, 1e-9);   // Cl
    EXPECT_NEAR(pt.GetElectronegativity(53), 2.66, 1e-9);   // I
    EXPECT_NEAR(pt.GetElectronegativity(55), 0.79, 1e-9);   // Cs (least electronegative tabulated)
    EXPECT_NEAR(pt.GetElectronegativity(82), 2.33, 1e-9);   // Pb
    EXPECT_EQ  (pt.GetElectronegativity(10), 0.0);          // Ne (noble gas: no ion)
    // F is the most electronegative; Cs/Fr-region the least -- the ordering the IonicSAD heuristic relies on.
    EXPECT_GT(pt.GetElectronegativity(9), pt.GetElectronegativity(11));   // F > Na  => F-, Na+
    EXPECT_GT(pt.GetElectronegativity(53), pt.GetElectronegativity(55));  // I > Cs  => I-, Cs+
}

TEST_F(PeriodicTableTests,ReadSaito)
{
    std::filesystem::path data_dir = COMMON_DATA_PATH;
    std::ifstream file(data_dir / "saito.json");
    assert(file);
    nlohmann::json jsondata;
    file >> jsondata;
    std::vector<ElementRecordSaito> pt;
    for (auto e:jsondata)
        pt.push_back(ElementRecordSaito(e));
    for (auto e:pt)
        std::cout << e;

}