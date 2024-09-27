// File: AtomTests.C  Test the DFT calculation for a hydrogen atom

#include "DFTTester.H"
#include "HartreeFockTester.H"
#include "DFTDataBase/DFTDataBase.H"
#include "Cluster/Atom.H"
#include "Misc/PeriodicTable.H"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

DFTDataBase  theDataBase("Atom.db");
PeriodicTable thePeriodicTable;
double eps_ro=1e-5; //Converge criterial for delta ro (charge density)
double eps_e=1e-4;

TEST_P(HartreeFockAtomTester, AtomsHFPolarized)
{
    int Z=GetParam();
    std::cout << "Testing atom " << thePeriodicTable.GetSymbol(Z) << std::endl;
    Init(new Atom(Z,0,Vector3D<double>(0,0,0)),thePeriodicTable.GetMaxL(Z),thePeriodicTable.GetNumUnpairedElectrons(Z));
    Iterate(1.0,eps_ro,40);
    double E_HF=thePeriodicTable.GetEnergyHF(Z);
    double error=fabs((E_HF-TotalEnergy())/E_HF);
    std::cout.precision(5);
    std::cout << "E_HF relative error=" << error*100.0 << "%" << std::endl;
    EXPECT_LT(error,eps_e);
};

TEST_P(DFTAtomTester, AtomsDFTPolarized)
{
    int Z=GetParam();
    std::cout << "Testing atom " << thePeriodicTable.GetSymbol(Z) << std::endl;
    Init(new Atom(Z,0,Vector3D<double>(0,0,0)),thePeriodicTable.GetSlaterAlpha(Z),thePeriodicTable.GetMaxL(Z),thePeriodicTable.GetNumUnpairedElectrons(Z));
    Iterate(1.0,eps_ro,40);
    double E_DFT=thePeriodicTable.GetEnergyDFT(Z);
    double error=fabs((E_DFT-TotalEnergy())/E_DFT);
    std::cout.precision(5);
    std::cout << "E_DFT relative error=" << error*100.0 << "%" << std::endl;
    EXPECT_LT(error,eps_e);
};

TEST_P(SemiHartreeFockAtomTester, AtomsSemiDFTPolarized)
{
    int Z=GetParam();
    std::cout << "Testing atom " << thePeriodicTable.GetSymbol(Z) << std::endl;
    Init(new Atom(Z,0,Vector3D<double>(0,0,0)),thePeriodicTable.GetSlaterAlpha(Z),thePeriodicTable.GetMaxL(Z),thePeriodicTable.GetNumUnpairedElectrons(Z));
    Iterate(1.0,eps_ro,20);
    double E_DFT=thePeriodicTable.GetEnergyDFT(Z);
    double error=fabs((E_DFT-TotalEnergy())/E_DFT);
    std::cout.precision(5);
    std::cout << "E_DFT relative error=" << error*100.0 << "%" << std::endl;
    EXPECT_LT(error,eps_e);
};

INSTANTIATE_TEST_CASE_P(AtomsHFPolarized,
                        HartreeFockAtomTester,
                        ::testing::Values(1,4,5,7,10,12,25));

INSTANTIATE_TEST_CASE_P(AtomsDFTPolarized,
                        DFTAtomTester,
                        ::testing::Values(1,4,5,7,10,12,25));

INSTANTIATE_TEST_CASE_P(AtomsSemiDFTPolarized,
                        SemiHartreeFockAtomTester,
                        ::testing::Range(3,17));

//TEST_F(HartreeFockAtomTester, UraniumPolarized)
//{
//    Init(new Atom(thePeriodicTable.GetZ("U"),0,Vector3D<double>(0,0,0)),3,4.0);
//    Iterate(1.0,eps_ro,50);
//    double expected_energy=-25658.417889;
//    EXPECT_LT(fabs((expected_energy-TotalEnergy())/expected_energy),eps_e);
//}


//TEST_F(DFTAtomTester, PeriodicTable)
//{
//    EXPECT_EQ(thePeriodicTable.GetSymbol(11),std::string("Na"));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(12),std::string("Mg"));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(13),std::string("Al"));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(14),std::string("Si"));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(15),std::string("P "));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(16),std::string("S "));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(17),std::string("Cl"));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(18),std::string("Ar"));
//    EXPECT_EQ(thePeriodicTable.GetSymbol(35),std::string("Br"));
//    EXPECT_EQ(thePeriodicTable.GetZ("Na"),11);
//    EXPECT_EQ(thePeriodicTable.GetZ("Mg"),12);
//    EXPECT_EQ(thePeriodicTable.GetZ("Al"),13);
//    EXPECT_EQ(thePeriodicTable.GetZ("Si"),14);
//    EXPECT_EQ(thePeriodicTable.GetZ("P" ),15);
//    EXPECT_EQ(thePeriodicTable.GetZ("S" ),16);
//    EXPECT_EQ(thePeriodicTable.GetZ("Cl"),17);
//    EXPECT_EQ(thePeriodicTable.GetZ("Ar"),18);
//    EXPECT_EQ(thePeriodicTable.GetZ("Br"),35);
//}
