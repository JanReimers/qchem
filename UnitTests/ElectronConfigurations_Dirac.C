// File: UnitTests/ElectronConfigurations_Dirac.C  Test relativistic atom electron configurations.


#include "gtest/gtest.h"
#include "tabulate/table.hpp"
#include <nlohmann/json.hpp>
#include <iostream>

import Common.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
import qchem.Symmetry.Atom_Dirac_EC;
// import qchem.Factory;
import qchem.Common.Strings;

using std::cout;
using std::endl;

// struct ml_Breakdown
// {
//     std::vector<int> ml_paired;     //List of ml values for paired orbitals
//     std::vector<int> ml_unpaired;   //List of ml values for unpaired orbitals
//     std::vector<int> ml_unoccupied; //List of ml values for empty orbitals
// };

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class Dirac_EC_Tests : public ::testing::Test
{
public:
    typedef sym_t sym_t;
    Dirac_EC_Tests() {}
    sym_t qn(int κ                   ) const {return Symmetry::ΩFactory(κ    );}
    sym_t qn(int κ, const rvec_t& mjs) const {return Symmetry::ΩFactory(κ,mjs);}
    
    static int GetN(const Atom_Dirac_EC& ac, const sym_t& sym, Spin s) 
    {
        Irrep qns(s,sym);
        return ac.GetN(qns);
    }
    // ml_Breakdown GetBreadown(const Atom_Dirac_EC& ec,size_t l) const;
    PeriodicTableSaito pt;
};

// ml_Breakdown Dirac_EC_Tests::GetBreadown(const Atom_Dirac_EC& ec,size_t l) const
// {
//     ec.itsNs.DebugCheck(); //Check self consistency
//     ml_Breakdown mls;
//     size_t g=2*l+1; //degenracy
//     size_t Nunp=abs(ec.itsNs.Nu[l]); //These can be negative Z=58 Ce.
//     size_t Npairs= g-Nunp; // For full shell systems this ends up being zero.
//     size_t Nempty=g-Npairs-Nunp;
//     if (Nempty==g) //Fix up for full shell systems.
//     {
//         Npairs=g;
//         Nempty=0;
//     }
//     // if (itsNs.Nu[l]>0)
//     {
//         int ml=-(int)l;
//         for (size_t i=0;i<Nunp  ;i++) mls.ml_unpaired  .push_back(ml++);
//         for (size_t i=0;i<Npairs;i++) mls.ml_paired    .push_back(ml++);
//         for (size_t i=0;i<Nempty;i++) mls.ml_unoccupied.push_back(ml++);
//         assert(ml==(int)(l+1));

//     }
    
//     return mls;
// }


// TEST_F(Dirac_EC_Tests, ElectronConfigurations)
// {
//     tabulate::Table BS_table;
//     BS_table.format().multi_byte_characters(true);
//     BS_table.add_row({"Z","Name","maxL","Nunpaired","Elconfig","s","p","d","f"});
//     for (size_t Z=1;Z<=92;Z++)
//     {
//         // cout << "Z=" << Z << endl;
//         Atom_Dirac_EC ec(Z);
//         tabulate::RowStream rs;
//         rs << Z;
//         rs << pt.GetSymbol(Z);
//         rs << pt.GetMaxL(Z);
//         rs << pt.GetNumUnpairedElectrons(Z);
//         const size_t* v=pt.GetValanceConfiguration(Z);
//         std::ostringstream os;
//         for (size_t i=0;i<4;i++)
//         {
//             if (v[i]>0)
//                 os << SPDFG[i] << superscripts[v[i]];
//         }
//         rs << os.str();       
//         for (size_t l=0;l<=ec.GetLMax();l++)
//         {
//             std::ostringstream os1;
//             ml_Breakdown ml=GetBreadown(ec,l);
//             size_t nunp=ml.ml_paired.size();
//             if (nunp!=0 && nunp!=2*l+1)
//                 for (size_t m=0;m<nunp;m++) os1 << "↑↓";
//             if (nunp>0) os1 << " ";
//             for (size_t m=0;m<ml.ml_unpaired.size();m++) os1 << "↑";
//             os1 << std::ends;
//             rs << os1.str();
//         }
//         BS_table.add_row(rs);
//     }
//     BS_table.column(5).format().font_color(l_colors[0]);
//     BS_table.column(6).format().font_color(l_colors[1]);
//     BS_table.column(7).format().font_color(l_colors[2]);
//     BS_table.column(8).format().font_color(l_colors[3]);
//     cout << BS_table << endl;  
// }

// TEST_F(Dirac_EC_Tests, BasisSets)
// {
//     tabulate::Table BS_table;
//     BS_table.format().multi_byte_characters(true);
//     BS_table.add_row({"Z","Name","s","p","d","f"});
//     nlohmann::json js = {{"N", 10},{"emin", 0.1},{"emax", 5000.0}
//                         ,{"type",abs_t::Slater}};
//     for (size_t Z=1;Z<=92;Z++)
//     {
//         tabulate::RowStream rs;
//         rs << Z;
//         rs << pt.GetSymbol(Z);
//         Real_BS* bs=BasisSet::Atom::Factory(js,Z);
//         size_t l=0;
//         std::ostringstream os[4];
//         for (auto ibs:bs->Iterate<Real_OIBS>())
//         {
//             size_t l1=Getl(ibs->GetSymmetry());
            // if (l>0 && l1==l) os[l] << endl;
            // if (l1>l) os[l++] << std::ends;
            // os[l] << sym;
//         }
//         for (auto& osl:os) rs << osl.str();
//         BS_table.add_row(rs);
//     }
//     BS_table.column(2).format().font_color(l_colors[0]);
//     BS_table.column(3).format().font_color(l_colors[1]);
//     BS_table.column(4).format().font_color(l_colors[2]);
//     BS_table.column(5).format().font_color(l_colors[3]);
//     cout << BS_table << endl;  
// }

TEST_F(Dirac_EC_Tests, Hydrogen)
{
    Atom_Dirac_EC ec(1);
    sym_t s12=qn(-1); //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),0);
}
TEST_F(Dirac_EC_Tests, Helium)
{
    Atom_Dirac_EC ec(2);
    sym_t s12=qn(-1); //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),1);
}
TEST_F(Dirac_EC_Tests, Lithium)
{
    Atom_Dirac_EC ec(3);
    sym_t s12=qn(-1);  //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),1);
}
TEST_F(Dirac_EC_Tests, Beryllium)
{
    Atom_Dirac_EC ec(4);
    sym_t s12=qn(-1);  //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
}
TEST_F(Dirac_EC_Tests, Boron)
{
    Atom_Dirac_EC ec(5);
    sym_t s12=qn(-1);//s+
    sym_t p12=qn(1); //p- j=1/2 this level is already s12lit.
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),0);
}
TEST_F(Dirac_EC_Tests, Carbon)
{
    Atom_Dirac_EC ec(6);
    sym_t s12=qn(-1);//s+
    sym_t p12=qn( 1); //p- //no need for mj s12litting
    sym_t p32=qn(-2,{-1.5}); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),0);
}
TEST_F(Dirac_EC_Tests, Nitrogen)
{
    Atom_Dirac_EC ec(7);
    sym_t s12=qn(-1);//s+
    sym_t p12=qn( 1);  //p- or p1/2  half full no mj splitting.
    sym_t p32=qn(-2); //p+ or p3/2 half full no mj splitting.
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),0);
}
TEST_F(Dirac_EC_Tests, Oxygen)
{
    Atom_Dirac_EC ec(8);
    sym_t s12=qn(-1);//s+
    sym_t p12=qn( 1); //p- or p1/2 full so no mj splitting.
    sym_t p32=qn(-2); //p+ or p3/2 half full so no mj splitting.
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),0);
}
TEST_F(Dirac_EC_Tests, Flourine)
{
    Atom_Dirac_EC ec(9);
    sym_t s12 =qn(-1); //s+
    sym_t p12 =qn( 1); //p- 
    sym_t p32p=qn(-2,{-1.5,-0.5}); //p+ or p3/2 paired
    sym_t p32u=qn(-2,{0.5}); //p+ or p3/2 on unpaired
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32p)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32p)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32u)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32u)),0);
}
TEST_F(Dirac_EC_Tests, Neon)
{
    Atom_Dirac_EC ec(10);
    sym_t s12=qn(-1); //s+
    sym_t p12=qn( 1); //p-
    sym_t p32=qn(-2); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
}
TEST_F(Dirac_EC_Tests, Sodium)
{
    Atom_Dirac_EC ec(11);
    sym_t s12=qn(-1);
    sym_t p12=qn( 1); //p-
    sym_t p32=qn(-2); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
}
TEST_F(Dirac_EC_Tests, Magnesium)
{
    Atom_Dirac_EC ec(12);
    sym_t s12=qn(-1);
    sym_t p12=qn( 1); //p-
    sym_t p32=qn(-2); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
}
TEST_F(Dirac_EC_Tests, Aluminium)
{
    Atom_Dirac_EC ec(13);
    sym_t s12=qn(-1);
    sym_t p12=qn( 1); //p-
    sym_t p32=qn(-2); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
}
// TEST_F(Dirac_EC_Tests, Silicon)
// {
//     Atom_Dirac_EC ec(14);
//     sym_t s12=qn(-1); //s+
//     sym_t p12=qn( 1); //p- //no need for mj s12litting
//     sym_t p32p=qn(-2,{-1.5,-0.5,0.5}); //p+
//     sym_t p32u=qn(-2,{ 1.5}); //p+
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
// }
// TEST_F(Dirac_EC_Tests, Phosphorus)
// {
//     Atom_Dirac_EC ec(15);
//     sym_t s12=qn(-1); //s+
//     sym_t p12=qn( 1); //p- or p1/2 half full no mj splitting.
//     sym_t p32=qn(-2); //p+ or p3/2 half full no mj splitting.
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),3);
// }
// TEST_F(Dirac_EC_Tests, Sulpher)
// {
//     Atom_Dirac_EC ec(16);
//     sym_t s12=qn(-1);//s+
//     sym_t p12=qn( 1); //p- or p1/2 full so no mj splitting.
//     sym_t p32=qn(-2); //p+ or p3/2 half full so no mj splitting.
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),3);
// }
// TEST_F(Dirac_EC_Tests, Chlorine)
// {
//     Atom_Dirac_EC ec(17);
//     sym_t s12 =qn(-1); //s+
//     sym_t p12 =qn( 1); //p- 
//     sym_t p32p=qn(-2,{-1.5,-0.5}); //p+ or p3/2 paired
//     sym_t p32u=qn(-2,{0.5}); //p+ or p3/2 on unpaired
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32p)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32p)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32u)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32u)),3);
// }
// TEST_F(Dirac_EC_Tests, Argon)
// {
//     Atom_Dirac_EC ec(18);
//     sym_t s12=qn(-1); //s+
//     sym_t p12=qn( 1); //p-
//     sym_t p32=qn(-2); //p+
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),4);
// }
// TEST_F(Dirac_EC_Tests, Potassium)
// {
//     Atom_Dirac_EC ec(19);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
// }
// TEST_F(Dirac_EC_Tests, Calcium)
// {
//     Atom_Dirac_EC ec(20);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
// }
// TEST_F(Dirac_EC_Tests, Scandium)
// {
//     Atom_Dirac_EC ec(21);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκmj(2,{-2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),1);
// }
// TEST_F(Dirac_EC_Tests, Titanium)
// {
//     Atom_Dirac_EC ec(22);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκmj(2,{-2,-1}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),2);
// }
// TEST_F(Dirac_EC_Tests, Vanadium)
// {
//     Atom_Dirac_EC ec(23);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκmj(2,{-2,-1,0}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),3);
// }
// TEST_F(Dirac_EC_Tests, Chromium)
// {
//     Atom_Dirac_EC ec(24);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Manganese)
// {
//     Atom_Dirac_EC ec(25); //| s²d⁵     |   |        | ↑↑↑↑↑     | 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Iron)
// {
//     Atom_Dirac_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1,0,1}));
//     sym_t d2(new Ωκmj(2,{2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
// }
// TEST_F(Dirac_EC_Tests, Cobalt)
// {
//     Atom_Dirac_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1,0}));
//     sym_t d2(new Ωκmj(2,{1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
// }
// TEST_F(Dirac_EC_Tests, Nickel)
// {
//     Atom_Dirac_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1}));
//     sym_t d2(new Ωκmj(2,{0,1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
// }
// TEST_F(Dirac_EC_Tests, Cop32er)
// {
//     Atom_Dirac_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Zinc)
// {
//     Atom_Dirac_EC ec(30); // | s²d¹⁰    |   |        |           |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Galium)
// {
//     Atom_Dirac_EC ec(31); // | s²p¹d¹⁰  |   | ↑      |           | 
//     sym_t s=qn(-1);
//     sym_t p1(new Ωκmj(1,{-1}));
//     sym_t p2(new Ωκmj(1,{0,1}));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Germanium)
// {
//     Atom_Dirac_EC ec(32); // | s²p²d¹⁰  |   | ↑↑     |
//     sym_t s=qn(-1);
//     sym_t p1(new Ωκmj(1,{-1,0}));
//     sym_t p2(new Ωκmj(1,{1}));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Arsenic)
// {
//     Atom_Dirac_EC ec(33); // | s²p³d¹⁰  |   | ↑↑↑    |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Selenium)
// {
//     Atom_Dirac_EC ec(34); //| s²p⁴d¹⁰  |   | ↑↓ ↑↑  |
//     sym_t s=qn(-1);
//     sym_t p_paired  (new Ωκmj(1,{1}));
//     sym_t p_unpaired(new Ωκmj(1,{-1,0}));
//     sym_t d(new Ωκ(2));

//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p_paired)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p_paired)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p_unpaired)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p_unpaired)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Bromine)
// {
//     Atom_Dirac_EC ec(35); // | s²p⁵d¹⁰  |   | ↑↓↑↓ ↑ | 
//     sym_t s=qn(-1);
//     sym_t p_paired  (new Ωκmj(1,{0,1}));
//     sym_t p_unpaired(new Ωκmj(1,{-1}));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p_paired)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p_paired)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p_unpaired)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p_unpaired)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Krypton)
// {
//     Atom_Dirac_EC ec(36);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Rubidium)
// {
//     Atom_Dirac_EC ec(37); //| s¹       | ↑ | 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Strontium)
// {
//     Atom_Dirac_EC ec(38); //| s²       |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Yttrium)
// {
//     Atom_Dirac_EC ec(39); // s²d¹     |   |        | ↑   
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2}));
//     sym_t d2(new Ωκmj(2,{-1,0,1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),1+1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),1+0);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),4);
// }
// TEST_F(Dirac_EC_Tests, Zirconium)
// {
//     Atom_Dirac_EC ec(40); //| s²d²     |   |        | ↑↑ 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1}));
//     sym_t d2(new Ωκmj(2,{0,1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),2+2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),2+0);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),3);
// }
// TEST_F(Dirac_EC_Tests, Niobium)
// {
//     Atom_Dirac_EC ec(41); //| s¹d⁴     | ↑ |        | ↑↑↑↑ 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1,0,1}));
//     sym_t d2(new Ωκmj(2,{2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),4+4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),4+0);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),1);
// }
// TEST_F(Dirac_EC_Tests, Molybdenum)
// {
//     Atom_Dirac_EC ec(42); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5+5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }
// TEST_F(Dirac_EC_Tests, Technesium)
// {
//     Atom_Dirac_EC ec(43); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5+5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }

// TEST_F(Dirac_EC_Tests, Iron)
// {
//     Atom_Dirac_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1,0,1}));
//     sym_t d2(new Ωκmj(2,{2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
// }

// TEST_F(Dirac_EC_Tests, Cobalt)
// {
//     Atom_Dirac_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1,0}));
//     sym_t d2(new Ωκmj(2,{1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
// }

// TEST_F(Dirac_EC_Tests, Nickel)
// {
//     Atom_Dirac_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d1(new Ωκmj(2,{-2,-1}));
//     sym_t d2(new Ωκmj(2,{0,1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),2);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
// }

// TEST_F(Dirac_EC_Tests, Cop32er)
// {
//     Atom_Dirac_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
//     sym_t d(new Ωκ(2));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
// }

