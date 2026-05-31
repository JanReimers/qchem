// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "tabulate/table.hpp"
#include <nlohmann/json.hpp>
#include <iostream>

import Common.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Ylm;
import qchem.Symmetry.Yl;
import qchem.Symmetry.AtomEC;
import qchem.Factory;
import qchem.Common.Strings;

using std::cout;
using std::endl;

struct ml_Breakdown
{
    std::vector<int> ml_paired;     //List of ml values for paired orbitals
    std::vector<int> ml_unpaired;   //List of ml values for unpaired orbitals
    std::vector<int> ml_unoccupied; //List of ml values for empty orbitals
};

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
    
    static int GetN(const Atom_EC& ac, const sym_t& sym, Spin s) 
    {
        Irrep_QNs qns(s,sym);
        return ac.GetN(qns);
    }
    ml_Breakdown GetBreadown(const Atom_EC& ec,size_t l) const;
    PeriodicTableSaito pt;
};

ml_Breakdown ElectronConfigurationTests::GetBreadown(const Atom_EC& ec,size_t l) const
{
    ec.itsNs.DebugCheck(); //Check self consistency
    ml_Breakdown mls;
    size_t g=2*l+1; //degenracy
    size_t Nunp=abs(ec.itsNs.Nu[l]); //These can be negative Z=58 Ce.
    size_t Npairs= g-Nunp; // For full shell systems this ends up being zero.
    size_t Nempty=g-Npairs-Nunp;
    if (Nempty==g) //Fix up for full shell systems.
    {
        Npairs=g;
        Nempty=0;
    }
    // if (itsNs.Nu[l]>0)
    {
        int ml=-(int)l;
        for (size_t i=0;i<Nunp  ;i++) mls.ml_unpaired  .push_back(ml++);
        for (size_t i=0;i<Npairs;i++) mls.ml_paired    .push_back(ml++);
        for (size_t i=0;i<Nempty;i++) mls.ml_unoccupied.push_back(ml++);
        assert(ml==(int)(l+1));

    }
    
    return mls;
}


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
        const size_t* v=pt.GetValanceConfiguration(Z);
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
            ml_Breakdown ml=GetBreadown(ec,l);
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
    BS_table.column(5).format().font_color(l_colors[0]);
    BS_table.column(6).format().font_color(l_colors[1]);
    BS_table.column(7).format().font_color(l_colors[2]);
    BS_table.column(8).format().font_color(l_colors[3]);
    cout << BS_table << endl;  

}

TEST_F(ElectronConfigurationTests, BasisSets)
{
    tabulate::Table BS_table;
    BS_table.format().multi_byte_characters(true);
    BS_table.add_row({"Z","Name","s","p","d","f"});
    
    nlohmann::json js = {{"N", 10},{"emin", 0.1},{"emax", 5000.0}
                        ,{"type",abs_t::Slater}};
    for (size_t Z=1;Z<=92;Z++)
    {
        tabulate::RowStream rs;
        rs << Z;
        rs << pt.GetSymbol(Z);
        Real_BS* bs=BasisSet::Atom::Factory(js,Z);
        size_t l=0;
        std::ostringstream os[4];
        for (auto ibs:bs->Iterate<Real_OIBS>())
        {
            const Angular_Sym& sym=ibs->CastSymmetry<Angular_Sym>();
            if (l>0 && sym.GetL()==l) os[l] << endl;
            if (sym.GetL()>l) os[l++] << std::ends;
            os[l] << sym;
        }
        for (auto& osl:os) rs << osl.str();
        BS_table.add_row(rs);
    }
    BS_table.column(2).format().font_color(l_colors[0]);
    BS_table.column(3).format().font_color(l_colors[1]);
    BS_table.column(4).format().font_color(l_colors[2]);
    BS_table.column(5).format().font_color(l_colors[3]);
    cout << BS_table << endl;  
}

TEST_F(ElectronConfigurationTests, Hydrogen)
{
    Atom_EC ec(1);
    sym_t s(new Yl_Sym(0));
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),0);
}
TEST_F(ElectronConfigurationTests, Helium)
{
    Atom_EC ec(2);
    sym_t s(new Yl_Sym(0));
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),1);
}
TEST_F(ElectronConfigurationTests, Lithium)
{
    Atom_EC ec(3);
    sym_t s(new Yl_Sym(0));
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),1);
}
TEST_F(ElectronConfigurationTests, Beryllium)
{
    Atom_EC ec(4);
    sym_t s(new Yl_Sym(0));
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
}
TEST_F(ElectronConfigurationTests, Boron)
{
    Atom_EC ec(5);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Ylm_Sym(1,{-1}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),0);
}
TEST_F(ElectronConfigurationTests, Carbon)
{
    Atom_EC ec(6);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Ylm_Sym(1,{-1,0}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),0);
}
TEST_F(ElectronConfigurationTests, Nitrogen)
{
    Atom_EC ec(7);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),0);
}
TEST_F(ElectronConfigurationTests, Oxygen)
{
    Atom_EC ec(8);
    sym_t s(new Yl_Sym(0));
    sym_t p_paired  (new Ylm_Sym(1,{1}));
    sym_t p_unpaired(new Ylm_Sym(1,{-1,0}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_paired)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_paired)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_unpaired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_unpaired)),0);
}
TEST_F(ElectronConfigurationTests, Flourine)
{
    Atom_EC ec(9);
    sym_t s(new Yl_Sym(0));
    sym_t p_paired  (new Ylm_Sym(1,{0,1}));
    sym_t p_unpaired(new Ylm_Sym(1,{-1}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_paired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_paired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_unpaired)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_unpaired)),0);
}
TEST_F(ElectronConfigurationTests, Neon)
{
    Atom_EC ec(10);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Sodium)
{
    Atom_EC ec(11);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Magnesium)
{
    Atom_EC ec(12);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Aluminium)
{
    Atom_EC ec(13);
    sym_t s(new Yl_Sym(0));
    sym_t p1(new Ylm_Sym(1,{-1}));
    sym_t p2(new Ylm_Sym(1,{0,1}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p1)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p2)),2);
}
TEST_F(ElectronConfigurationTests, Silicon)
{
    Atom_EC ec(14);
    sym_t s(new Yl_Sym(0));
    sym_t p1(new Ylm_Sym(1,{-1,0}));
    sym_t p2(new Ylm_Sym(1,{1}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p1)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p2)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p2)),1);
}
TEST_F(ElectronConfigurationTests, Phosphorus)
{
    Atom_EC ec(15);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Sulpher)
{
    Atom_EC ec(16);
    sym_t s(new Yl_Sym(0));
    sym_t p_paired  (new Ylm_Sym(1,{1}));
    sym_t p_unpaired(new Ylm_Sym(1,{-1,0}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_paired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_paired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_unpaired)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_unpaired)),2);
}
TEST_F(ElectronConfigurationTests, Chlorine)
{
    Atom_EC ec(17);
    sym_t s(new Yl_Sym(0));
    sym_t p_paired  (new Ylm_Sym(1,{0,1}));
    sym_t p_unpaired(new Ylm_Sym(1,{-1}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_paired)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_paired)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_unpaired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_unpaired)),1);
}
TEST_F(ElectronConfigurationTests, Argon)
{
    Atom_EC ec(18);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
}
TEST_F(ElectronConfigurationTests, Potassium)
{
    Atom_EC ec(19);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
}
TEST_F(ElectronConfigurationTests, Calcium)
{
    Atom_EC ec(20);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
}
TEST_F(ElectronConfigurationTests, Scandium)
{
    Atom_EC ec(21);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Ylm_Sym(2,{-2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),1);
}
TEST_F(ElectronConfigurationTests, Titanium)
{
    Atom_EC ec(22);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Ylm_Sym(2,{-2,-1}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),2);
}
TEST_F(ElectronConfigurationTests, Vanadium)
{
    Atom_EC ec(23);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Ylm_Sym(2,{-2,-1,0}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),3);
}
TEST_F(ElectronConfigurationTests, Chromium)
{
    Atom_EC ec(24);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
}
TEST_F(ElectronConfigurationTests, Manganese)
{
    Atom_EC ec(25); //| s²d⁵     |   |        | ↑↑↑↑↑     | 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
}
TEST_F(ElectronConfigurationTests, Iron)
{
    Atom_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d1(new Ylm_Sym(2,{-2,-1,0,1}));
    sym_t d2(new Ylm_Sym(2,{2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),0);
}
TEST_F(ElectronConfigurationTests, Cobalt)
{
    Atom_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d1(new Ylm_Sym(2,{-2,-1,0}));
    sym_t d2(new Ylm_Sym(2,{1,2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),0);
}
TEST_F(ElectronConfigurationTests, Nickel)
{
    Atom_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d1(new Ylm_Sym(2,{-2,-1}));
    sym_t d2(new Ylm_Sym(2,{0,1,2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),0);
}
TEST_F(ElectronConfigurationTests, Copper)
{
    Atom_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Zinc)
{
    Atom_EC ec(30); // | s²d¹⁰    |   |        |           |
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Galium)
{
    Atom_EC ec(31); // | s²p¹d¹⁰  |   | ↑      |           | 
    sym_t s(new Yl_Sym(0));
    sym_t p1(new Ylm_Sym(1,{-1}));
    sym_t p2(new Ylm_Sym(1,{0,1}));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p1)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p2)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p2)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Germanium)
{
    Atom_EC ec(32); // | s²p²d¹⁰  |   | ↑↑     |
    sym_t s(new Yl_Sym(0));
    sym_t p1(new Ylm_Sym(1,{-1,0}));
    sym_t p2(new Ylm_Sym(1,{1}));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p1)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p1)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Arsenic)
{
    Atom_EC ec(33); // | s²p³d¹⁰  |   | ↑↑↑    |
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Selenium)
{
    Atom_EC ec(34); //| s²p⁴d¹⁰  |   | ↑↓ ↑↑  |
    sym_t s(new Yl_Sym(0));
    sym_t p_paired  (new Ylm_Sym(1,{1}));
    sym_t p_unpaired(new Ylm_Sym(1,{-1,0}));
    sym_t d(new Yl_Sym(2));

    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_paired)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_paired)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_unpaired)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_unpaired)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Bromine)
{
    Atom_EC ec(35); // | s²p⁵d¹⁰  |   | ↑↓↑↓ ↑ | 
    sym_t s(new Yl_Sym(0));
    sym_t p_paired  (new Ylm_Sym(1,{0,1}));
    sym_t p_unpaired(new Ylm_Sym(1,{-1}));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_paired)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_paired)),6);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p_unpaired)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p_unpaired)),2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Krypton)
{
    Atom_EC ec(36);
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Rubidium)
{
    Atom_EC ec(37); //| s¹       | ↑ | 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Strontium)
{
    Atom_EC ec(38); //| s²       |
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Yttrium)
{
    Atom_EC ec(39); // s²d¹     |   |        | ↑   
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d1(new Ylm_Sym(2,{-2}));
    sym_t d2(new Ylm_Sym(2,{-1,0,1,2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),1+1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),1+0);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),4);
}
TEST_F(ElectronConfigurationTests, Zirconium)
{
    Atom_EC ec(40); //| s²d²     |   |        | ↑↑ 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d1(new Ylm_Sym(2,{-2,-1}));
    sym_t d2(new Ylm_Sym(2,{0,1,2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),2+2);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),2+0);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),3);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),3);
}
TEST_F(ElectronConfigurationTests, Niobium)
{
    Atom_EC ec(41); //| s¹d⁴     | ↑ |        | ↑↑↑↑ 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d1(new Ylm_Sym(2,{-2,-1,0,1}));
    sym_t d2(new Ylm_Sym(2,{2}));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),4+4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),4+0);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),1);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),1);
}
TEST_F(ElectronConfigurationTests, Molybdenum)
{
    Atom_EC ec(42); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5+5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Technesium)
{
    Atom_EC ec(43); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
    sym_t s(new Yl_Sym(0));
    sym_t p(new Yl_Sym(1));
    sym_t d(new Yl_Sym(2));
    
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5+5);
    EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
}

// TEST_F(ElectronConfigurationTests, Iron)
// {
//     Atom_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
//     sym_t s(new Yl_Sym(0));
//     sym_t p(new Yl_Sym(1));
//     sym_t d1(new Ylm_Sym(2,{-2,-1,0,1}));
//     sym_t d2(new Ylm_Sym(2,{2}));
    
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),1);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),0);
// }

// TEST_F(ElectronConfigurationTests, Cobalt)
// {
//     Atom_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
//     sym_t s(new Yl_Sym(0));
//     sym_t p(new Yl_Sym(1));
//     sym_t d1(new Ylm_Sym(2,{-2,-1,0}));
//     sym_t d2(new Ylm_Sym(2,{1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),2);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),2);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),3);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),0);
// }

// TEST_F(ElectronConfigurationTests, Nickel)
// {
//     Atom_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
//     sym_t s(new Yl_Sym(0));
//     sym_t p(new Yl_Sym(1));
//     sym_t d1(new Ylm_Sym(2,{-2,-1}));
//     sym_t d2(new Ylm_Sym(2,{0,1,2}));
    
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d2)),3);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d1)),2);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d1)),0);
// }

// TEST_F(ElectronConfigurationTests, Copper)
// {
//     Atom_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
//     sym_t s(new Yl_Sym(0));
//     sym_t p(new Yl_Sym(1));
//     sym_t d(new Yl_Sym(2));
    
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,s)),4);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,s)),3);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,p)),6);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Up  ,d)),5);
//     EXPECT_EQ(ec.GetN(Irrep_QNs(Spin::Down,d)),5);
// }

