// File: UnitTests/ElectronConfigurations.C  Test non-relativistic atom electron configurations.


#include "gtest/gtest.h"
#include <stdexcept>
#include "tabulate/table.hpp"
#include <iostream>

import qchem.PeriodicTable;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.Factory;

import qchem.ElectronConfiguration.AtomNR;
import qchem.ElectronConfiguration.Molecule;   // Molecule_EC (spin-native molecular occupation, OpenWork B3)
import qchem.Strings;
using namespace qchem;

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
    typedef sym_t sym_t;
    ElectronConfigurationTests() {}
    sym_t qn(int l) const {return Symmetry::YFactory(l);}
    sym_t qn(int l, const ivec_t& mls) const {return Symmetry::YFactory(l,mls);}
    
    static int GetN(const Atom_EC& ac, const sym_t& sym, Spin s) 
    {
        Irrep qns(s,sym);
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

// TEST_F(ElectronConfigurationTests, BasisSets)
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
//         Real_BS* bs=BasisSet::Radial::Factory(js,Z);
//         size_t l=0;
//         std::ostringstream os[4];
//         for (auto ibs:bs->Iterate<Real_OIBS>())
//         {
//             size_t l1=Getl(ibs->GetSymmetry());
//             if (l>0 && l1==l) os[l] << endl;
//             if (l1>l) os[l++] << std::ends;
//             os[l] << ibs->GetSymmetry();
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

TEST_F(ElectronConfigurationTests, Hydrogen)
{
    Atom_EC ec(1);
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),0);
}
TEST_F(ElectronConfigurationTests, Helium)
{
    Atom_EC ec(2);
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),1);
}
TEST_F(ElectronConfigurationTests, Lithium)
{
    Atom_EC ec(3);
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),1);
}
TEST_F(ElectronConfigurationTests, Beryllium)
{
    Atom_EC ec(4);
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
}
TEST_F(ElectronConfigurationTests, Boron)
{
    Atom_EC ec(5);
    sym_t s=qn(0);
    sym_t p=qn(1,{-1});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),0);
}
TEST_F(ElectronConfigurationTests, Carbon)
{
    Atom_EC ec(6);
    sym_t s=qn(0);
    sym_t p=qn(1,{-1,0});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),0);
}
TEST_F(ElectronConfigurationTests, Nitrogen)
{
    Atom_EC ec(7);
    sym_t s=qn(0);
    sym_t p=qn(1);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),0);
}
TEST_F(ElectronConfigurationTests, Oxygen)
{
    Atom_EC ec(8);
    sym_t s=qn(0);
    sym_t p1  =qn(1,{1}); //paired
    sym_t p2=qn(1,{-1,0}); //unpaired
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),0);
}
TEST_F(ElectronConfigurationTests, Flourine)
{
    Atom_EC ec(9);
    sym_t s=qn(0);
    sym_t p1  =qn(1,{0,1}); //paired
    sym_t p2=qn(1,{-1}); //unpaired
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),0);
}
TEST_F(ElectronConfigurationTests, Neon)
{
    Atom_EC ec(10);
    sym_t s=qn(0);
    sym_t p=qn(1);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Sodium)
{
    Atom_EC ec(11);
    sym_t s=qn(0);
    sym_t p=qn(1);

    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Magnesium)
{
    Atom_EC ec(12);
    sym_t s=qn(0);
    sym_t p=qn(1);

    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Aluminium)
{
    Atom_EC ec(13);
    sym_t s=qn(0);
    sym_t p1=qn(1,{-1});
    sym_t p2=qn(1,{0,1});

    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),2);
}
TEST_F(ElectronConfigurationTests, Silicon)
{
    Atom_EC ec(14);
    sym_t s=qn(0);
    sym_t p1=qn(1,{-1,0});
    sym_t p2=qn(1,{1});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),1);
}
// Pseudo-atom (valence-only) config: Si keeps just its 3s^2 3p^2 valence over the s and p atomic irreps,
// with NO core and NO folding of the filled 3s^2 into the core (contrast the all-electron Silicon test
// above, where the s irrep carries 1s+2s+3s = 6 electrons).  Also exercises the Spin::None path for an
// EMPTY (virtual) irrep -- the p{+1} sub-irrep of the open p^2 shell -- which must return 0, not garbage.
TEST_F(ElectronConfigurationTests, SiliconPseudoValence)
{
    PseudoAtom_EC ec(14);
    sym_t s     = qn(0);
    sym_t p_occ = qn(1,{-1,0});   // the two singly-occupied p orbitals (Hund's rule, spin up)
    sym_t p_emp = qn(1,{1});      // the empty p sub-irrep

    // Spin-polarized occupations: s^2 (1 up, 1 down) + p^2 (2 up in {-1,0}, 0 down).
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s    )),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s    )),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p_occ)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p_occ)),0);

    // Un-polarized totals (Spin::None): 2 in s, 2 in p{-1,0}, and 0 in the empty p{+1} irrep.
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,s    )),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,p_occ)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,p_emp)),0);   // empty-irrep lookup must be 0 (was UB)
}
// Pseudo-ion configs: the net charge adjusts the neutral pseudo-valence.  F- (charge -1) adds one electron
// to close the p shell (2s^2 2p^5 -> 2s^2 2p^6); Na+ (charge +1) removes the lone 3s electron (empty valence).
TEST_F(ElectronConfigurationTests, FluorideAnionPseudoValence)
{
    PseudoAtom_EC ec(9, -1);          // F-
    sym_t s=qn(0), p=qn(1);           // closed p shell -> single (un-split) p irrep
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)), 3);   // p^6 closed: 3 up + 3 down
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)), 3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,s)), 2);   // 8 valence electrons total (2s^2 2p^6)
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,p)), 6);
}
TEST_F(ElectronConfigurationTests, SodiumCationPseudoValence)
{
    PseudoAtom_EC ec(11, +1);         // Na+: removes the 3s^1 -> empty valence
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,s)), 0);   // 0 valence electrons (closed Ne-like core in the PP)
}
// TRANSITION-METAL cation: the physical (aufbau-reverse) removal order empties the HIGHER-n 4s shell
// before 3d, so Mn2+ (q7 valence 4s^2 3d^5 - 2) is 3d^5 S=5/2 -- the MnO/LiMn2O4 local moment the
// spin-SAD seed stores (doc/SCFSeedingPlan.md sec 10).  NOT 4s^2 3d^3 (the old highest-l-first rule).
TEST_F(ElectronConfigurationTests, ManganeseCationPseudoValence)
{
    PseudoAtom_EC ec(25, +2);         // Mn2+: removes 4s^2 -> d^5 half-filled, 5 unpaired
    sym_t s=qn(0), d=qn(2);           // half-filled d: single (un-split) d irrep
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)), 5);   // Hund d^5: all majority
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)), 0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,d)), 5);
}
// And the multi-open-shell TM cation: Mn+ = 4s^1 3d^5 (the ^7S ground state) -- the per-shell Hund sum
// counts BOTH open shells (1+5 = 6 unpaired).
TEST_F(ElectronConfigurationTests, ManganeseMonocationPseudoValence)
{
    PseudoAtom_EC ec(25, +1);         // Mn+: removes ONE 4s electron
    sym_t s=qn(0), d=qn(2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)), 5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)), 0);
}
// SEMICORE pseudopotentials (2026-10-31): a q above the neutral valence keeps closed core shells -- the outermost of the noble-gas
// core, highest n first.  Until then q was ignored and a q9 Na atom carried ONE electron (the periodic table's neutral valence).
TEST_F(ElectronConfigurationTests, SodiumSemicoreQ9)
{
    sym_t s=qn(0), p=qn(1);
    PseudoAtom_EC na(11, 0, 9);       // Na q9: 2s^2 2p^6 3s^1
    EXPECT_EQ(na.GetN(Irrep(Spin::None,s)), 3);
    EXPECT_EQ(na.GetN(Irrep(Spin::None,p)), 6);
    PseudoAtom_EC nap(11, +1, 9);     // Na+ q9: the cation empties 3s, never the semicore -> 2s^2 2p^6
    EXPECT_EQ(nap.GetN(Irrep(Spin::None,s)), 2);
    EXPECT_EQ(nap.GetN(Irrep(Spin::None,p)), 6);
    EXPECT_EQ(nap.GetN(Irrep(Spin::Up,p)), 3);
    EXPECT_EQ(nap.GetN(Irrep(Spin::Down,p)), 3);
}
TEST_F(ElectronConfigurationTests, LithiumSemicoreQ3)
{
    PseudoAtom_EC li(3, 0, 3);        // Li q3: 1s^2 2s^1 -- both in l=0
    EXPECT_EQ(li.GetN(Irrep(Spin::None,qn(0))), 3);
    PseudoAtom_EC lip(3, +1, 3);      // Li+ q3: 1s^2
    EXPECT_EQ(lip.GetN(Irrep(Spin::None,qn(0))), 2);
}
TEST_F(ElectronConfigurationTests, ManganeseSemicoreQ15)
{
    PseudoAtom_EC mn(25, +2, 15);     // Mn2+ q15: 3s^2 3p^6 3d^5, the d^5 still all majority
    sym_t s=qn(0), p=qn(1), d=qn(2);
    EXPECT_EQ(mn.GetN(Irrep(Spin::None,s)), 2);
    EXPECT_EQ(mn.GetN(Irrep(Spin::None,p)), 6);
    EXPECT_EQ(mn.GetN(Irrep(Spin::Up  ,d)), 5);
    EXPECT_EQ(mn.GetN(Irrep(Spin::Down,d)), 0);
}
TEST_F(ElectronConfigurationTests, ASemicoreQThatSplitsAShellIsRefused)
{
    EXPECT_THROW(PseudoAtom_EC(11, 0, 5), std::invalid_argument);   // 4 semicore electrons: half of 2p
    EXPECT_THROW(PseudoAtom_EC(11, 0, 0), std::invalid_argument);   // fewer than the neutral valence
    EXPECT_NO_THROW(PseudoAtom_EC(11, 0, 1));                       // the valence-only q, explicitly
}
TEST_F(ElectronConfigurationTests, Phosphorus)
{
    Atom_EC ec(15);
    sym_t s=qn(0);
    sym_t p=qn(1);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),3);
}
TEST_F(ElectronConfigurationTests, Sulpher)
{
    Atom_EC ec(16);
    sym_t s=qn(0);
    sym_t p1=qn(1,{1}); //paired
    sym_t p2=qn(1,{-1,0});  //unpaired
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),2);
}
TEST_F(ElectronConfigurationTests, Chlorine)
{
    Atom_EC ec(17);
    sym_t s=qn(0);
    sym_t p1=qn(1,{0,1}); //paired
    sym_t p2=qn(1,{-1}); //unpaired
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),1);
}
TEST_F(ElectronConfigurationTests, Argon)
{
    Atom_EC ec(18);
    sym_t s=qn(0);
    sym_t p=qn(1);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
}
TEST_F(ElectronConfigurationTests, Potassium)
{
    Atom_EC ec(19);
    sym_t s=qn(0);
    sym_t p=qn(1);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
}
TEST_F(ElectronConfigurationTests, Calcium)
{
    Atom_EC ec(20);
    sym_t s=qn(0);
    sym_t p=qn(1);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
}
TEST_F(ElectronConfigurationTests, Scandium)
{
    Atom_EC ec(21);
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2,{-2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),1);
}
TEST_F(ElectronConfigurationTests, Titanium)
{
    Atom_EC ec(22);
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2,{-2,-1});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),2);
}
TEST_F(ElectronConfigurationTests, Vanadium)
{
    Atom_EC ec(23);
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2,{-2,-1,0});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),3);
}
TEST_F(ElectronConfigurationTests, Chromium)
{
    Atom_EC ec(24);
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
}
TEST_F(ElectronConfigurationTests, Manganese)
{
    Atom_EC ec(25); //| s²d⁵     |   |        | ↑↑↑↑↑     | 
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
}
TEST_F(ElectronConfigurationTests, Iron)
{
    Atom_EC ec(26); //| s²d⁶     |   |        | ↑↓ ↑↑↑↑   |
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d1=qn(2,{-2,-1,0,1});
    sym_t d2=qn(2,{2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
}
TEST_F(ElectronConfigurationTests, Cobalt)
{
    Atom_EC ec(27); // | s²d⁷     |   |        | ↑↓↑↓ ↑↑↑  | 
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d1=qn(2,{-2,-1,0});
    sym_t d2=qn(2,{1,2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
}
TEST_F(ElectronConfigurationTests, Nickel)
{
    Atom_EC ec(28); // | s²d⁸     |   |        | ↑↓↑↓↑↓ ↑↑ |
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d1=qn(2,{-2,-1});
    sym_t d2=qn(2,{0,1,2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),0);
}
TEST_F(ElectronConfigurationTests, Copper)
{
    Atom_EC ec(29); // | s¹d¹⁰    | ↑ |        |           |
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Zinc)
{
    Atom_EC ec(30); // | s²d¹⁰    |   |        |           |
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Galium)
{
    Atom_EC ec(31); // | s²p¹d¹⁰  |   | ↑      |           | 
    sym_t s=qn(0);
    sym_t p1=qn(1,{-1});
    sym_t p2=qn(1,{0,1});
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Germanium)
{
    Atom_EC ec(32); // | s²p²d¹⁰  |   | ↑↑     |
    sym_t s=qn(0);
    sym_t p1=qn(1,{-1,0});
    sym_t p2=qn(1,{1});
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Arsenic)
{
    Atom_EC ec(33); // | s²p³d¹⁰  |   | ↑↑↑    |
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Selenium)
{
    Atom_EC ec(34); //| s²p⁴d¹⁰  |   | ↑↓ ↑↑  |
    sym_t s=qn(0);
    sym_t p1  =qn(1,{1});
    sym_t p2=qn(1,{-1,0});
    sym_t d=qn(2);

    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Bromine)
{
    Atom_EC ec(35); // | s²p⁵d¹⁰  |   | ↑↓↑↓ ↑ | 
    sym_t s=qn(0);
    sym_t p1  =qn(1,{0,1});
    sym_t p2=qn(1,{-1});
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p1)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p1)),6);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p2)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p2)),2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Krypton)
{
    Atom_EC ec(36);
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Rubidium)
{
    Atom_EC ec(37); //| s¹       | ↑ | 
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Strontium)
{
    Atom_EC ec(38); //| s²       |
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Yttrium)
{
    Atom_EC ec(39); // s²d¹     |   |        | ↑   
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d1=qn(2,{-2});
    sym_t d2=qn(2,{-1,0,1,2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),1+1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),1+0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),4);
}
TEST_F(ElectronConfigurationTests, Zirconium)
{
    Atom_EC ec(40); //| s²d²     |   |        | ↑↑ 
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d1=qn(2,{-2,-1});
    sym_t d2=qn(2,{0,1,2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),2+2);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),2+0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),3);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),3);
}
TEST_F(ElectronConfigurationTests, Niobium)
{
    Atom_EC ec(41); //| s¹d⁴     | ↑ |        | ↑↑↑↑ 
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d1=qn(2,{-2,-1,0,1});
    sym_t d2=qn(2,{2});
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d1)),4+4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d1)),4+0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d2)),1);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d2)),1);
}
TEST_F(ElectronConfigurationTests, Molybdenum)
{
    Atom_EC ec(42); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑ 
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);
    
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5+5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}
TEST_F(ElectronConfigurationTests, Technesium)
{
    Atom_EC ec(43); //| s¹d⁵     | ↑ |        | ↑↑↑↑↑
    sym_t s=qn(0);
    sym_t p=qn(1);
    sym_t d=qn(2);

    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)),5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,p)),9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,d)),5+5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,d)),5);
}

// --- Molecule_EC: spin-native molecular occupation (OpenWork B3) ----------------------------------
// The single-arg Molecule_EC(Ne) is the minimal-spin collapse (singlet for even Ne, doublet for odd) and
// must reproduce the historical closed-shell behaviour byte-for-byte; the two-arg form is the spin-native
// open-shell (nUp, nDown).  Molecule_EC ignores the spatial irrep (global aufbau fills per spin channel),
// so any sym exercises GetN's per-spin dispatch.
TEST_F(ElectronConfigurationTests, MoleculeClosedShellEven)
{
    Molecule_EC ec(10);                 // 10 e- -> singlet: 5 up, 5 down (== the old closed-shell split)
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,s)),10);
}
TEST_F(ElectronConfigurationTests, MoleculeMinimalSpinOdd)
{
    Molecule_EC ec(9);                  // 9 e- -> doublet: 5 up, 4 down (the historical odd-Ne collapse)
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 5);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 4);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,s)), 9);
}
TEST_F(ElectronConfigurationTests, MoleculeOpenShellTriplet)
{
    Molecule_EC ec(9, 7);               // explicit open shell: 9 up, 7 down (triplet, nUp-nDown = 2)
    sym_t s=qn(0);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Up  ,s)), 9);
    EXPECT_EQ(ec.GetN(Irrep(Spin::Down,s)), 7);
    EXPECT_EQ(ec.GetN(Irrep(Spin::None,s)),16);
}


//----------------------------------------------------------------------------------------------
//  A STATED configuration (2026-11-15, OpenWork row "All-electron ION seeds").
namespace
{
int Spin2S(const Atom_EC& ec)   // 2S from the occupations: up minus down over every irrep
{
    int up=0, dn=0;
    for (auto ir : ec.GetIrreps()) { up+=ec.GetN(Irrep(Spin::Up,ir)); dn+=ec.GetN(Irrep(Spin::Down,ir)); }
    return up-dn;
}
int Total(const Atom_EC& ec)
{
    int n=0;
    for (auto ir : ec.GetIrreps()) n+=ec.GetN(Irrep(Spin::Up,ir))+ec.GetN(Irrep(Spin::Down,ir));
    return n;
}
}

TEST(StatedConfiguration, NISTsNotationParsesWithItsCoreExpanded)
{
    const auto fe=ParseConfiguration("[Ar] 3d^6 4s^1");
    ASSERT_TRUE(fe) << fe.Error();
    EXPECT_EQ(fe.Value().coreSymbol, "Ar");
    EXPECT_EQ(fe.Value().CoreElectrons(), 18);
    EXPECT_EQ(fe.Value().NamedElectrons(), 7);
    EXPECT_EQ(ConfigurationText(fe.Value().All()), "1s^2 2s^2 2p^6 3s^2 3p^6 3d^6 4s^1");
    const auto bare=ParseConfiguration("1s2 2s2 2p6");   // the ^ is optional, no core
    ASSERT_TRUE(bare) << bare.Error();
    EXPECT_EQ(bare.Value().CoreElectrons(), 0);
    EXPECT_EQ(bare.Value().NamedElectrons(), 10);
    for (const char* bad : {"[Xy] 3d^6", "[Ar 3d^6", "3x^2", "3d^", "2d^1", "3d^11", "[Ar] 3p^2", "4s^1 4s^1", ""})
        EXPECT_FALSE(ParseConfiguration(bad)) << "\"" << bad << "\" must be refused";
}

// The radial atom fills each l block from its lowest level up, so only the highest occupied n of each l may be open and no level may be
// skipped; an all-electron atom also starts every l at n = l+1.
TEST(StatedConfiguration, TheRadialAtomHoldsOnlyBottomUpFilledLBlocks)
{
    auto shells=[](const char* t){ return ParseConfiguration(t).Value().All(); };
    EXPECT_EQ(WhyNotRepresentable(shells("[Ar] 3d^6 4s^1"), true), "");
    EXPECT_NE(WhyNotRepresentable(shells("[Ne] 3s^1 5s^1"), true), "") << "a hole at 4s";
    EXPECT_NE(WhyNotRepresentable(shells("1s^2 2s^1 3s^1"), true), "") << "2s open below an occupied 3s";
    EXPECT_NE(WhyNotRepresentable(shells("2s^2"), true), "") << "an all-electron atom starts s at 1s";
    EXPECT_EQ(WhyNotRepresentable(shells("3d^5"), false), "") << "a pseudo-atom's d block starts at its valence";
}

// Fe+ is NIST's [Ar] 3d^6 4s^1 (high spin: four unpaired d and the s), where its isoelectronic neutral Mn is 3d^5 4s^2 (2S = 5 too, but
// another configuration).  The stated EC carries 25 electrons, 2S = 5, and reads back the text it was stated as.
TEST(StatedConfiguration, IronCationIsItsOwnConfigurationNotManganeses)
{
    const Atom_EC fe(ParseConfiguration("[Ar] 3d^6 4s^1").Value().All());
    EXPECT_EQ(Total(fe), 25);
    EXPECT_EQ(Spin2S(fe), 5);
    EXPECT_EQ(ConfigurationText(fe.Shells()), "1s^2 2s^2 2p^6 3s^2 3p^6 3d^6 4s^1");
    const Atom_EC mn(25);
    EXPECT_EQ(ConfigurationText(mn.Shells()), "1s^2 2s^2 2p^6 3s^2 3p^6 3d^5 4s^2") << "the derived rule: the isoelectronic neutral";
    const Atom_EC v(ParseConfiguration("[Ar] 3d^4").Value().All());   // V+: 4s moved INTO 3d
    EXPECT_EQ(Spin2S(v), 4);
    EXPECT_THROW(Atom_EC(ParseConfiguration("[Ne] 3s^1 5s^1").Value().All()), std::logic_error) << "an unchecked shape is a broken invariant";
}

// A PSEUDO-ATOM's configuration text names its PHYSICAL shells, over whatever core its pseudopotential replaces -- not only a Madelung
// prefix of Z - q electrons.  Kr q8 replaces [Ar] 3d^10, Cs q9 [Kr] 4d^10: Shells() used to fill the core in aufbau order and THREW for
// both (and for every lanthanide, actinide and 4d semicore q), so the atom facade could not run them (2026-11-15, row 189).
TEST(StatedConfiguration, APseudoAtomNamesItsShellsOverTheCoreItsPseudopotentialReplaces)
{
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(14, 0, 4).Shells()),  "3s^2 3p^2");
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(25, 0, 15).Shells()), "3s^2 3p^6 3d^5 4s^2") << "Mn q15: the semicore 3s 3p";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(25, 3, 7).Shells()),  "3d^4")               << "Mn3+ q7: the cation empties 4s first";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(36, 0, 8).Shells()),  "4s^2 4p^6")          << "Kr q8 over [Ar] 3d^10";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(55, 0, 9).Shells()),  "5s^2 5p^6 6s^1")     << "Cs q9 over [Kr] 4d^10";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(38, 0, 10).Shells()), "4s^2 4p^6 5s^2")     << "Sr q10 over [Ar] 3d^10";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(55, 0, 1).Shells()),  "6s^1")               << "Cs q1 over [Xe]";
}

// RULE B (controller ruling 2026-11-15, row 189): a pseudopotential whose q is SMALLER than the periodic table's valence holds that valence's
// closed f^14, then d^10, subshells in its core.  Before it, every such q was refused ("keeps FEWER electrons"), so the default smallest-q
// set had no pseudo-atom for I, Ga, Zn, Cu or any 5d metal.  A q no whole-shell fold reaches is refused, never guessed.
TEST(StatedConfiguration, ASmallQPseudopotentialFoldsTheClosedSubshellsItReplacesIntoItsCore)
{
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(53, 0, 7).Shells()),  "5s^2 5p^5")          << "I q7 over [Kr] 4d^10";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(31, 0, 3).Shells()),  "4s^2 4p^1")          << "Ga q3 over [Ar] 3d^10";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(31, 0, 13).Shells()), "3d^10 4s^2 4p^1")    << "Ga q13 keeps its 3d";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(30, 0, 2).Shells()),  "4s^2")               << "Zn q2";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(29, 0, 1).Shells()),  "4s^1")               << "Cu q1";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(74, 0, 6).Shells()),  "5d^4 6s^2")          << "W q6 over [Xe] 4f^14";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(79, 0, 1).Shells()),  "6s^1")               << "Au q1 folds 4f^14 and 5d^10";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(72, 0, 12).Shells()), "5s^2 5p^6 5d^2 6s^2") << "Hf q12 folds 4f^14, keeps the 5s 5p semicore";
    EXPECT_EQ(ConfigurationText(PseudoAtom_EC(53, -1, 7).Shells()), "5s^2 5p^6")          << "I- q7";
    EXPECT_THROW(PseudoAtom_EC(26, 0, 4), std::invalid_argument) << "Fe q4: 3d^6 4s^2 has no closed subshell to fold";
}

// THE TWO ROUTES AGREE: every neutral atom's derived EC, read back as subshells and stated again, occupies every irrep identically
// (the stated route is the derived one wherever the derived one is right).
TEST_F(ElectronConfigurationTests, EveryNeutralStatedAgainOccupiesEveryIrrepAsTheDerivedOne)
{
    for (int Z=1; Z<=92; Z++)
    {
        const Atom_EC derived(Z);
        const ShellOccupations shells=derived.Shells();
        ASSERT_EQ(WhyNotRepresentable(shells, true), "") << "Z=" << Z << " " << ConfigurationText(shells);
        const Atom_EC stated(shells);
        EXPECT_EQ(Total(stated), Z) << "Z=" << Z;
        EXPECT_EQ(Spin2S(stated), Spin2S(derived)) << "Z=" << Z << " " << ConfigurationText(shells);
        for (auto ir : derived.GetIrreps())
            for (Spin s : {Spin::Up, Spin::Down})
                EXPECT_EQ(stated.GetN(Irrep(s,ir)), derived.GetN(Irrep(s,ir))) << "Z=" << Z << " " << ConfigurationText(shells);
    }
}
