// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/Misc/PeriodicTable.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include <Spin.H>
#include <iostream>

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class ElectronConfigurationTests : public ::testing::Test
{
public:
    ElectronConfigurationTests() {}
};

TEST_F(ElectronConfigurationTests, Ntotal)
{
//    int Z=44;
    for (int Z=1;Z<=94;Z++)
    {
        AtomElectronConfiguration ec(Z);
        EXPECT_EQ(ec.GetN(),Z);
    }
}

TEST_F(ElectronConfigurationTests, SumSpin)
{
//    int Z=89;
    for (int Z=1;Z<=94;Z++)
    {
        AtomElectronConfiguration ec(Z);
        cout << "Z=" << Z << endl;
        EXPECT_EQ(ec.GetN(),ec.GetN(Spin::Up)+ec.GetN(Spin::Down));
    }
}

TEST_F(ElectronConfigurationTests, SumL)
{
//    int Z=89;
    for (int Z=1;Z<=94;Z++)
    {
        AtomElectronConfiguration ec(Z);
        cout << "Z=" << Z << endl;
        int Nl=0;
        for (int l=0;l<=3;l++)
        {
            SphericalSymmetryQN qn(l);
            Nl+=ec.GetN(qn);
        }
        EXPECT_EQ(ec.GetN(),Nl);
    }
}

TEST_F(ElectronConfigurationTests, SumLAndSpin)
{
//    int Z=31;
    for (int Z=1;Z<=56;Z++)
    {
        AtomElectronConfiguration ec(Z);
        cout << "Z=" << Z << endl;
        int Nlu=0;
        int Nld=0;
        for (int l=0;l<=3;l++)
        {
            SphericalSymmetryQN qn(l);
            int nlu=ec.GetN(qn,Spin::Up);
            int nld=ec.GetN(qn,Spin::Down);
            cout << "l,nu,nd = " << l << " " << nlu << " " << nld << endl;
            EXPECT_EQ(ec.GetN(qn),nlu+nld);
            Nlu+=nlu;
            Nld+=nld;
        }
        EXPECT_EQ(ec.GetN(Spin::Up),Nlu);
        EXPECT_EQ(ec.GetN(Spin::Down),Nld);
        EXPECT_EQ(ec.GetN(),Nlu+Nld);
    }
}
