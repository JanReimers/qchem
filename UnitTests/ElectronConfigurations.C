// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/BasisSet/Atom/EC.H"
#include "Imp/Misc/PeriodicTable.H"
#include "Imp/BasisSet/Atom/ml/Ylm.H"
#include <Irrep_QNs.H>
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
    typedef Irrep_QNs::sym_t sym_t;
    ElectronConfigurationTests() {}
    sym_t qn(int l) const {return sym_t(new Yl_Sym(l));}
    sym_t qn(int l, const std::vector<int>& ml) const {return sym_t(new Ylm_Sym(l,ml));}
    
    static int GetN(const Atom_EC& ac) {return ac.GetN();}
    static int GetN(const Atom_EC& ac, Spin s) {return ac.GetN(s);}
    static int GetN(const Atom_EC& ac, const sym_t& s) {return ac.GetN(s);}
    static int GetN(const Atom_EC& ac, const sym_t& sym, Spin s) 
    {
        Irrep_QNs qns(s,sym);
        return ac.GetN(qns);
    }
};

TEST_F(ElectronConfigurationTests, Ntotal)
{
//    int Z=44;
    for (int Z=1;Z<=94;Z++)
    {
        Atom_EC ec(Z);
        EXPECT_EQ(GetN(ec),Z);
    }
}

TEST_F(ElectronConfigurationTests, SumSpin)
{
//    int Z=89;
    for (int Z=1;Z<=94;Z++)
    {
        Atom_EC ec(Z);
//        cout << "Z=" << Z << endl;
        EXPECT_EQ(GetN(ec),GetN(ec,Spin::Up)+GetN(ec,Spin::Down));
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
            Nl+=GetN(ec,qn(l));
        }
        EXPECT_EQ(GetN(ec),Nl);
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
            auto yl=qn(l);
            int nlu=GetN(ec,yl,Spin::Up);
            int nld=GetN(ec,yl,Spin::Down);
//            cout << "l,nu,nd = " << l << " " << nlu << " " << nld << endl;
            EXPECT_EQ(GetN(ec,yl),nlu+nld);
            Nlu+=nlu;
            Nld+=nld;
        }
        EXPECT_EQ(GetN(ec,Spin::Up),Nlu);
        EXPECT_EQ(GetN(ec,Spin::Down),Nld);
        EXPECT_EQ(GetN(ec),Nlu+Nld);
    }
}

TEST_F(ElectronConfigurationTests, SPconfigs)
{
    { //N
        Atom_EC ec(7);
        EXPECT_EQ(GetN(ec,Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,Spin::Down),2);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),2);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),2);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),0);
    }
    { //Flourine
        Atom_EC ec(9);
        EXPECT_EQ(GetN(ec,Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),2);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),2);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),2);
    }
    { //Sodium
        Atom_EC ec(11);
        EXPECT_EQ(GetN(ec,Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,Spin::Down),5);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),2);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),3);
    }
    { //Mg
        Atom_EC ec(12);
        EXPECT_EQ(GetN(ec,Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),3);
    }
    { //P
        Atom_EC ec(15);
        EXPECT_EQ(GetN(ec,Spin::Up  ),9);
        EXPECT_EQ(GetN(ec,Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),3);
    }
    { //As
        Atom_EC ec(33);
        EXPECT_EQ(GetN(ec,Spin::Up  ),18);
        EXPECT_EQ(GetN(ec,Spin::Down),15);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),9);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),5);
    }
    { //Sb
        Atom_EC ec(51);
        EXPECT_EQ(GetN(ec,Spin::Up  ),27);
        EXPECT_EQ(GetN(ec,Spin::Down),24);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),5);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),12);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),9);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),10);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),10);
    }
    { //Bi
        Atom_EC ec(83);
        EXPECT_EQ(GetN(ec,Spin::Up  ),43);
        EXPECT_EQ(GetN(ec,Spin::Down),40);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),12);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),15);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),7);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),7);
    }
}

TEST_F(ElectronConfigurationTests, Dconfigs)
{
    { //Sc
        Atom_EC ec(21);
        EXPECT_EQ(GetN(ec,Spin::Up  ),11);
        EXPECT_EQ(GetN(ec,Spin::Down),10);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),1);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Ti
        Atom_EC ec(22);
        EXPECT_EQ(GetN(ec,Spin::Up  ),12);
        EXPECT_EQ(GetN(ec,Spin::Down),10);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),2);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //V
        Atom_EC ec(23);
        EXPECT_EQ(GetN(ec,Spin::Up  ),13);
        EXPECT_EQ(GetN(ec,Spin::Down),10);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),3);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Cr
        Atom_EC ec(24);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),9);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Mn
        Atom_EC ec(25);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),10);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Fe
        Atom_EC ec(26);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),11);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),1);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Co
        Atom_EC ec(27);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),12);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),2);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Ni
        Atom_EC ec(28);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),13);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),3);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Cu
        Atom_EC ec(29);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),14);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),3);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),5);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }
    { //Zn
        Atom_EC ec(30);
        EXPECT_EQ(GetN(ec,Spin::Up  ),15);
        EXPECT_EQ(GetN(ec,Spin::Down),15);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Up  ),4);
        EXPECT_EQ(GetN(ec,qn(0),Spin::Down),4);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Up  ),6);
        EXPECT_EQ(GetN(ec,qn(1),Spin::Down),6);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Up  ),5);
        EXPECT_EQ(GetN(ec,qn(2),Spin::Down),5);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Up  ),0);
        EXPECT_EQ(GetN(ec,qn(3),Spin::Down),0);
    }

}

// TEST_F(ElectronConfigurationTests, Ylm_SP)
// {
//     for (int Z=4;Z<=10;Z++)
//     {
//         cout << "Z=" << Z << endl;
//         Atom_EC ec(Z);
//         for (int l=0;l<=3;l++)
//         {
// //            cout << "  l=" << l << endl;
//             int nlu=GetN(ec,qn(l),Spin::Up);
//             int nld=GetN(ec,qn(l),Spin::Down);
//             ml_Breakdown mls=ec.GetBreadown(l);
//             int nlu1=0,nld1=0;
//             nlu1+=GetN(ec,qn(l,mls.ml_unpaired),Spin::Up);
//             nlu1+=GetN(ec,qn(l,mls.ml_paired),Spin::Up);
//             nld1+=GetN(ec,qn(l,mls.ml_unpaired),Spin::Down);
//             nld1+=GetN(ec,qn(l,mls.ml_paired),Spin::Down);
//             EXPECT_EQ(nlu,nlu1);
//             EXPECT_EQ(nld,nld1);
//         }
//     }
    
// }

//TEST_F(ElectronConfigurationTests, Ylm_D)
//{
//    Atom_EC ec(41);
//    for (int l=0;l<=2;l++)
//    {
//        for (int m=-l;m<=l;m++)
//        {
//            cout << l << " " << m << " " << GetN(ec,qn(l,m),Spin::Up) << endl;
//            cout << l << " " << m << " " << GetN(ec,qn(l,m),Spin::Down) << endl;      
//        }        
//    }
//    
//}
