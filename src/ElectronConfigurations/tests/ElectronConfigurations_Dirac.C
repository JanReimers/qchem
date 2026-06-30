// File: UnitTests/ElectronConfigurations_Dirac.C  Test relativistic atom electron configurations.


#include "gtest/gtest.h"
#include "tabulate/table.hpp"
#include <iostream>

import qchem.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
import qchem.ElectronConfiguration.AtomDirac;
import qchem.Strings;
using namespace qchem;

using std::cout;
using std::endl;

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
    
    static int GetN(const AtomDirac_EC& ac, const sym_t& sym, Spin s) 
    {
        Irrep qns(s,sym);
        return ac.GetN(qns);
    }
    PeriodicTableSaito pt;
};


// TEST_F(Dirac_EC_Tests, ElectronConfigurations)
// {
//     tabulate::Table BS_table;
//     BS_table.format().multi_byte_characters(true);
//     BS_table.add_row({"Z","Name","maxL","Nunpaired","Elconfig","s","p","d","f"});
//     for (size_t Z=1;Z<=92;Z++)
//     {
//         // cout << "Z=" << Z << endl;
//         AtomDirac_EC ec(Z);
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

TEST_F(Dirac_EC_Tests, Hydrogen)
{
    AtomDirac_EC ec(1);
    sym_t s12=qn(-1); //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),0);
}
TEST_F(Dirac_EC_Tests, Helium)
{
    AtomDirac_EC ec(2);
    sym_t s12=qn(-1); //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),1);
}
TEST_F(Dirac_EC_Tests, Lithium)
{
    AtomDirac_EC ec(3);
    sym_t s12=qn(-1);  //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),1);
}
TEST_F(Dirac_EC_Tests, Beryllium)
{
    AtomDirac_EC ec(4);
    sym_t s12=qn(-1);  //s+
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
}
TEST_F(Dirac_EC_Tests, Boron)
{
    AtomDirac_EC ec(5);
    sym_t s12=qn(-1);//s+
    sym_t p12=qn(1); //p- j=1/2 this level is already s12lit.
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),0);
}
TEST_F(Dirac_EC_Tests, Carbon)
{
    AtomDirac_EC ec(6);
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
    AtomDirac_EC ec(7);
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
    AtomDirac_EC ec(8);
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
    AtomDirac_EC ec(9);
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
    AtomDirac_EC ec(10);
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
    AtomDirac_EC ec(11);
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
    AtomDirac_EC ec(12);
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
    AtomDirac_EC ec(13);
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
TEST_F(Dirac_EC_Tests, Silicon)
{
    AtomDirac_EC ec(14);
    sym_t s12=qn(-1); //s+
    sym_t p12=qn( 1); //p- //no need for mj s12litting
    sym_t p32p=qn(-2,{-1.5,-0.5,0.5}); //p+
    sym_t p32u=qn(-2,{ 1.5}); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32p)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32p)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32u)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32u)),0);
}
TEST_F(Dirac_EC_Tests, Phosphorus)
{
    AtomDirac_EC ec(15);
    sym_t s12=qn(-1); //s+
    sym_t p12=qn( 1); //p- or p1/2 half full no mj splitting.
    sym_t p32=qn(-2); //p+ or p3/2 half full no mj splitting.
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
}
TEST_F(Dirac_EC_Tests, Sulpher)
{
    AtomDirac_EC ec(16);
    sym_t s12=qn(-1);//s+
    sym_t p12=qn( 1); //p- or p1/2 full so no mj splitting.
    sym_t p32=qn(-2); //p+ or p3/2 half full so no mj splitting.
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),2);
}
TEST_F(Dirac_EC_Tests, Chlorine)
{
    AtomDirac_EC ec(17);
    sym_t s12 =qn(-1); //s+
    sym_t p12 =qn( 1); //p- 
    sym_t p32p=qn(-2,{-1.5,-0.5,0.5}); //p+ or p3/2 paired
    sym_t p32u=qn(-2,{1.5}); //p+ or p3/2 on unpaired
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32p)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32p)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32u)),0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32u)),1);
}
TEST_F(Dirac_EC_Tests, Argon)
{
    AtomDirac_EC ec(18);
    sym_t s12=qn(-1); //s+
    sym_t p12=qn( 1); //p-
    sym_t p32=qn(-2); //p+
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s12)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p12)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p32)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p32)),4);
}
// Closed-shell noble gases beyond Ar.  Each kappa-irrep accumulates all its n-levels,
// split half/half into Up/Down (Kramers).  GetN(Up,kappa)=GetN(Down,kappa)=N_levels*(2j+1)/2.
TEST_F(Dirac_EC_Tests, Xenon)
{
    AtomDirac_EC ec(54); // [Kr]4d10 5s2 5p6 ; 1s..5s(5) 2p..5p(4) 3d,4d(2)
    sym_t s12=qn(-1); //s1/2  5 levels
    sym_t p12=qn( 1); //p1/2  4 levels
    sym_t p32=qn(-2); //p3/2  4 levels
    sym_t d32=qn( 2); //d3/2  2 levels
    sym_t d52=qn(-3); //d5/2  2 levels

    EXPECT_EQ(GetN(ec,s12,Spin::Up),5); EXPECT_EQ(GetN(ec,s12,Spin::Down),5);
    EXPECT_EQ(GetN(ec,p12,Spin::Up),4); EXPECT_EQ(GetN(ec,p12,Spin::Down),4);
    EXPECT_EQ(GetN(ec,p32,Spin::Up),8); EXPECT_EQ(GetN(ec,p32,Spin::Down),8);
    EXPECT_EQ(GetN(ec,d32,Spin::Up),4); EXPECT_EQ(GetN(ec,d32,Spin::Down),4);
    EXPECT_EQ(GetN(ec,d52,Spin::Up),6); EXPECT_EQ(GetN(ec,d52,Spin::Down),6);
}
TEST_F(Dirac_EC_Tests, Radon)
{
    AtomDirac_EC ec(86); // [Xe]4f14 5d10 6s2 6p6 ; s(6) p(5) d(3) f(1)
    sym_t s12=qn(-1); //s1/2  6 levels
    sym_t p12=qn( 1); //p1/2  5 levels
    sym_t p32=qn(-2); //p3/2  5 levels
    sym_t d32=qn( 2); //d3/2  3 levels
    sym_t d52=qn(-3); //d5/2  3 levels
    sym_t f52=qn( 3); //f5/2  1 level
    sym_t f72=qn(-4); //f7/2  1 level

    EXPECT_EQ(GetN(ec,s12,Spin::Up),6);  EXPECT_EQ(GetN(ec,s12,Spin::Down),6);
    EXPECT_EQ(GetN(ec,p12,Spin::Up),5);  EXPECT_EQ(GetN(ec,p12,Spin::Down),5);
    EXPECT_EQ(GetN(ec,p32,Spin::Up),10); EXPECT_EQ(GetN(ec,p32,Spin::Down),10);
    EXPECT_EQ(GetN(ec,d32,Spin::Up),6);  EXPECT_EQ(GetN(ec,d32,Spin::Down),6);
    EXPECT_EQ(GetN(ec,d52,Spin::Up),9);  EXPECT_EQ(GetN(ec,d52,Spin::Down),9);
    EXPECT_EQ(GetN(ec,f52,Spin::Up),3);  EXPECT_EQ(GetN(ec,f52,Spin::Down),3);
    EXPECT_EQ(GetN(ec,f72,Spin::Up),4);  EXPECT_EQ(GetN(ec,f72,Spin::Down),4);
}
// TEST_F(Dirac_EC_Tests, Potassium)
// {
//     AtomDirac_EC ec(19);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
// }
// TEST_F(Dirac_EC_Tests, Calcium)
// {
//     AtomDirac_EC ec(20);
//     sym_t s=qn(-1);
//     sym_t p(new Ωκ(1));
    
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
// }
// TEST_F(Dirac_EC_Tests, Scandium)
// {
//     AtomDirac_EC ec(21);
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
//     AtomDirac_EC ec(22);
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
//     AtomDirac_EC ec(23);
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
//     AtomDirac_EC ec(24);
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
//     AtomDirac_EC ec(25); //| s²d⁵     |   |        | ↑↑↑↑↑     | 
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
//     AtomDirac_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
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
//     AtomDirac_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
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
//     AtomDirac_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
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
//     AtomDirac_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
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
//     AtomDirac_EC ec(30); // | s²d¹⁰    |   |        |           |
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
//     AtomDirac_EC ec(31); // | s²p¹d¹⁰  |   | ↑      |           | 
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
//     AtomDirac_EC ec(32); // | s²p²d¹⁰  |   | ↑↑     |
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
//     AtomDirac_EC ec(33); // | s²p³d¹⁰  |   | ↑↑↑    |
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
//     AtomDirac_EC ec(34); //| s²p⁴d¹⁰  |   | ↑↓ ↑↑  |
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
//     AtomDirac_EC ec(35); // | s²p⁵d¹⁰  |   | ↑↓↑↓ ↑ | 
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
//     AtomDirac_EC ec(36);
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
//     AtomDirac_EC ec(37); //| s¹       | ↑ | 
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
//     AtomDirac_EC ec(38); //| s²       |
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
//     AtomDirac_EC ec(39); // s²d¹     |   |        | ↑   
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
//     AtomDirac_EC ec(40); //| s²d²     |   |        | ↑↑ 
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
//     AtomDirac_EC ec(41); //| s¹d⁴     | ↑ |        | ↑↑↑↑ 
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
//     AtomDirac_EC ec(42); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
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
//     AtomDirac_EC ec(43); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
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
//     AtomDirac_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
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
//     AtomDirac_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
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
//     AtomDirac_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
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
//     AtomDirac_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
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

