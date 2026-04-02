// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "tabulate/table.hpp"
#include <nlohmann/json.hpp>
#include <iostream>

import Common.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Ylm;
import qchem.Symmetry.AtomEC;
import qchem.BasisSet;
import qchem.Factory;
import qchem.WaveFunction.Internal.CompositeWF;

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
    PeriodicTable pt;
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
std::string superscripts[]={"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹","¹⁰","¹¹","¹²","¹³","¹⁴","¹⁵","¹⁶","¹⁷","¹⁸"};
TEST_F(ElectronConfigurationTests, ElectronConfigurations)
{
    
    
    tabulate::Table BS_table;
    BS_table.format().multi_byte_characters(true);
    BS_table.add_row({"Z","Name","maxL","Nunpaired","Elconfig","s","p","d","f"});

    for (size_t Z=1;Z<=92;Z++)
    {
        // cout << "Z=" << Z << endl;
        Atom_EC ec(Z);
       
        tabulate::RowStream rs;
        rs << Z;
        rs << pt.GetSymbol(Z);
        rs << pt.GetMaxL(Z);
        rs << pt.GetNumUnpairedElectrons(Z);
        int* v=pt.GetValanceConfiguration(Z);
        std::ostringstream os;
        for (size_t i=0;i<4;i++)
        {
            if (v[i]>0)
                os << SPDFG[i] << superscripts[v[i]];
        }
        rs << os.str();

        
        for (size_t l=0;l<=ec.GetLMax();l++)
        {
            std::ostringstream os1;
            ml_Breakdown ml=ec.GetBreadown(l);
            size_t nunp=ml.ml_paired.size();
            if (nunp!=0 && nunp!=2*l+1)
                for (size_t m=0;m<nunp;m++) os1 << "↑↓";
            if (nunp>0) os1 << " ";
            for (size_t m=0;m<ml.ml_unpaired.size();m++) os1 << "↑";
            os1 << std::ends;
            rs << os1.str();
        }
        
        
        BS_table.add_row(rs);
    }
    BS_table.column(5).format().font_color(CompositeWF::l_colors[0]);
    BS_table.column(6).format().font_color(CompositeWF::l_colors[1]);
    BS_table.column(7).format().font_color(CompositeWF::l_colors[2]);
    BS_table.column(8).format().font_color(CompositeWF::l_colors[3]);
    cout << BS_table << endl;  

}

TEST_F(ElectronConfigurationTests, BasisSets)
{
    tabulate::Table BS_table;
    BS_table.format().multi_byte_characters(true);
    BS_table.add_row({"Z","Name","s","p","d","f"});
    
    nlohmann::json js = {{"N", 10},{"emin", 0.1},{"emax", 5000.0}};
    for (size_t Z=1;Z<=92;Z++)
    {
        tabulate::RowStream rs;
        rs << Z;
        rs << pt.GetSymbol(Z);
        BasisSet* bs=BasisSetAtom::Factory(BasisSetAtom::Type::Slater,js,Z);
        size_t l=0;
        std::ostringstream os[4];
        for (auto ibs:bs->Iterate<Real_OIBS>())
        {
            const Angular_Sym* sym=dynamic_cast<const Angular_Sym*>(ibs->GetSymmetry().get());
            if (l>0 && sym->GetL()==l) os[l] << endl;
            if (sym->GetL()>l) os[l++] << std::ends;
            os[l] << *sym;
        }
        for (auto& osl:os) rs << osl.str();
        BS_table.add_row(rs);
    }
    BS_table.column(2).format().font_color(CompositeWF::l_colors[0]);
    BS_table.column(3).format().font_color(CompositeWF::l_colors[1]);
    BS_table.column(4).format().font_color(CompositeWF::l_colors[2]);
    BS_table.column(5).format().font_color(CompositeWF::l_colors[3]);
    cout << BS_table << endl;  
}