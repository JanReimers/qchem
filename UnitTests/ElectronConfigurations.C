// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/Misc/PeriodicTable.H"
#include "Imp/Symmetry/YlmQN.H"
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
    YlQN qn(int l) const {return YlQN(l);}
    YlmQN qn(int l, int m) const {return YlmQN(l,m);}
};

TEST_F(ElectronConfigurationTests, Ntotal)
{
//    int Z=44;
    for (int Z=1;Z<=94;Z++)
    {
        Atom_EC ec(Z);
        EXPECT_EQ(ec.GetN(),Z);
    }
}

TEST_F(ElectronConfigurationTests, SumSpin)
{
//    int Z=89;
    for (int Z=1;Z<=94;Z++)
    {
        Atom_EC ec(Z);
//        cout << "Z=" << Z << endl;
        EXPECT_EQ(ec.GetN(),ec.GetN(Spin::Up)+ec.GetN(Spin::Down));
    }
}

TEST_F(ElectronConfigurationTests, SumL)
{
//    int Z=89;
    for (int Z=1;Z<=94;Z++)
    {
        Atom_EC ec(Z);
        // cout << "Z=" << Z << endl;
        int Nl=0;
        for (int l=0;l<=3;l++)
        {
            Nl+=ec.GetN(qn(l));
        }
        EXPECT_EQ(ec.GetN(),Nl);
    }
}

TEST_F(ElectronConfigurationTests, SumLAndSpin)
{
//    int Z=58;
    for (int Z=1;Z<=94;Z++)
    {
        Atom_EC ec(Z);
//        cout << "Z=" << Z << endl;
        int Nlu=0;
        int Nld=0;
        for (int l=0;l<=3;l++)
        {
            YlQN qn(l);
            int nlu=ec.GetN(qn,Spin::Up);
            int nld=ec.GetN(qn,Spin::Down);
//            cout << "l,nu,nd = " << l << " " << nlu << " " << nld << endl;
            EXPECT_EQ(ec.GetN(qn),nlu+nld);
            Nlu+=nlu;
            Nld+=nld;
        }
        EXPECT_EQ(ec.GetN(Spin::Up),Nlu);
        EXPECT_EQ(ec.GetN(Spin::Down),Nld);
        EXPECT_EQ(ec.GetN(),Nlu+Nld);
    }
}

TEST_F(ElectronConfigurationTests, SPconfigs)
{
    { //N
        Atom_EC ec(7);
        EXPECT_EQ(ec.GetN(Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(Spin::Down),2);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),2);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),2);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),0);
    }
    { //Flourine
        Atom_EC ec(9);
        EXPECT_EQ(ec.GetN(Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),2);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),2);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),2);
    }
    { //Sodium
        Atom_EC ec(11);
        EXPECT_EQ(ec.GetN(Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(Spin::Down),5);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),2);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),3);
    }
    { //Mg
        Atom_EC ec(12);
        EXPECT_EQ(ec.GetN(Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),3);
    }
    { //P
        Atom_EC ec(15);
        EXPECT_EQ(ec.GetN(Spin::Up  ),9);
        EXPECT_EQ(ec.GetN(Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),3);
    }
    { //As
        Atom_EC ec(33);
        EXPECT_EQ(ec.GetN(Spin::Up  ),18);
        EXPECT_EQ(ec.GetN(Spin::Down),15);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),9);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),5);
    }
    { //Sb
        Atom_EC ec(51);
        EXPECT_EQ(ec.GetN(Spin::Up  ),27);
        EXPECT_EQ(ec.GetN(Spin::Down),24);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),5);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),12);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),9);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),10);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),10);
    }
    { //Bi
        Atom_EC ec(83);
        EXPECT_EQ(ec.GetN(Spin::Up  ),43);
        EXPECT_EQ(ec.GetN(Spin::Down),40);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),12);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),15);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),7);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),7);
    }
}

TEST_F(ElectronConfigurationTests, Dconfigs)
{
    { //Sc
        Atom_EC ec(21);
        EXPECT_EQ(ec.GetN(Spin::Up  ),11);
        EXPECT_EQ(ec.GetN(Spin::Down),10);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),1);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Ti
        Atom_EC ec(22);
        EXPECT_EQ(ec.GetN(Spin::Up  ),12);
        EXPECT_EQ(ec.GetN(Spin::Down),10);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),2);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //V
        Atom_EC ec(23);
        EXPECT_EQ(ec.GetN(Spin::Up  ),13);
        EXPECT_EQ(ec.GetN(Spin::Down),10);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),3);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Cr
        Atom_EC ec(24);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),9);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Mn
        Atom_EC ec(25);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),10);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Fe
        Atom_EC ec(26);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),11);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),1);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Co
        Atom_EC ec(27);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),12);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),2);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Ni
        Atom_EC ec(28);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),13);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),3);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Cu
        Atom_EC ec(29);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),14);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),3);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),5);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }
    { //Zn
        Atom_EC ec(30);
        EXPECT_EQ(ec.GetN(Spin::Up  ),15);
        EXPECT_EQ(ec.GetN(Spin::Down),15);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Up  ),4);
        EXPECT_EQ(ec.GetN(qn(0),Spin::Down),4);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Up  ),6);
        EXPECT_EQ(ec.GetN(qn(1),Spin::Down),6);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Up  ),5);
        EXPECT_EQ(ec.GetN(qn(2),Spin::Down),5);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Up  ),0);
        EXPECT_EQ(ec.GetN(qn(3),Spin::Down),0);
    }

}

TEST_F(ElectronConfigurationTests, Ylm_SP)
{
    for (int Z=1;Z<=94;Z++)
    {
//        cout << "Z=" << Z << endl;
        Atom_EC ec(Z);
        for (int l=0;l<=3;l++)
        {
//            cout << "  l=" << l << endl;
            int nlu=ec.GetN(qn(l),Spin::Up);
            int nld=ec.GetN(qn(l),Spin::Down);
            int nlu1=0,nld1=0;
            for (int m=-l;m<=l;m++)
            {
                nlu1+=ec.GetN(qn(l,m),Spin::Up);
                nld1+=ec.GetN(qn(l,m),Spin::Down);
            }
            EXPECT_EQ(nlu,nlu1);
            EXPECT_EQ(nld,nld1);
        }
    }
    
}

//TEST_F(ElectronConfigurationTests, Ylm_D)
//{
//    Atom_EC ec(41);
//    for (int l=0;l<=2;l++)
//    {
//        for (int m=-l;m<=l;m++)
//        {
//            cout << l << " " << m << " " << ec.GetN(qn(l,m),Spin::Up) << endl;
//            cout << l << " " << m << " " << ec.GetN(qn(l,m),Spin::Down) << endl;      
//        }        
//    }
//    
//}
