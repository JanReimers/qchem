#include <memory>
#include <vector>
#include <ranges>
#include <algorithm>
#include <iterator>
#include "gtest/gtest.h"

using std::cout;
using std::endl;

import qchem.Structure;
import Structure.UnitCell;
import qchem.Lattice;
import qchem.Math;

//  The storage-agnostic index iterator must model a full forward iterator so
//  that range-based for AND the C++20 <ranges> adaptors both work on it.
static_assert(std::forward_iterator<Structure::const_iterator>);
static_assert(std::ranges::forward_range<Structure>);

class StructureTests : public ::testing::Test
{};

TEST_F(StructureTests, Atom)
{
    Atom B(5,0);
    Atom F_minus(9,-1);
    Atom Na_plus(11,1,{1,1,1});
    EXPECT_EQ(B.GetNumElectrons(),5);
    EXPECT_EQ(F_minus.GetNumElectrons(),10);
    EXPECT_EQ(Na_plus.GetNumElectrons(),10);
    cout << "           B=" << B << endl;
    cout << "          F-=" << F_minus << endl;
    cout << "         Na+=" << Na_plus << endl;
  
}

TEST_F(StructureTests, AtomValueSemantics)
{
    //  With the dummy/self-pointer kludge gone, an Atom is a plain value type:
    //  it can be stored by value and relocated (vector growth moves elements),
    //  and it still iterates as a one-element structure consisting of itself.
    Atom na(11,1,{1,2,3});
    std::vector<Atom> v;
    for (int i=0;i<8;i++) v.push_back(na); //Forces reallocation -> moves.
    v[2].itsR={7,8,9};
    EXPECT_EQ(v[0].itsR,rvec3_t(1,2,3));
    EXPECT_EQ(v[2].itsR,rvec3_t(7,8,9));
    EXPECT_EQ(v[2].GetNumAtoms(),1u);
    for (auto* a:v[2]) EXPECT_EQ(a->itsR,rvec3_t(7,8,9)); //self-range after moves.
}

TEST_F(StructureTests, RangeIteration)
{
    Molecule h2o;
    h2o.Insert(new Atom(8,0,{0,0,0}));
    h2o.Insert(new Atom(1,0,{0,0,1}));
    h2o.Insert(new Atom(1,0,{0,1,0}));

    //  Plain range-based for: yields Atom* by value through the index iterator.
    int n=0,Zsum=0;
    for (auto a:h2o) {n++; Zsum+=a->itsZ;}
    EXPECT_EQ(n,3);
    EXPECT_EQ(Zsum,10);

    //  C++20 <ranges> adaptors now compose on a Structure with no storage exposed.
    EXPECT_EQ(std::ranges::count_if(h2o,[](Atom* a){return a->itsZ==1;}),2);
    auto Zs=h2o | std::views::transform([](Atom* a){return a->itsZ;});
    EXPECT_EQ(std::ranges::max(Zs),8);
}

TEST_F(StructureTests, UnitCell)
{
    //  Geometry is compared with EXPECT_NEAR: A-matrix vs metric-tensor routes
    //  to the same value differ at the ULP level, which is physically irrelevant.
    UnitCell cubic(4,4,4,90,90,90);
    EXPECT_NEAR(cubic.GetCellVolume(),64,1e-9);
    EXPECT_NEAR(cubic.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(cubic.GetDistance(rvec3_t{1,1,1}),4*sqrt(3),1e-9);
    EXPECT_EQ(cubic.GetNumCells(9),Vector3D<int>(3,3,3));
    rvec3_t rc=cubic.ToCartesian({0.5,0.5,0.5}); // r = A f
    EXPECT_NEAR(rc.x,2,1e-12); EXPECT_NEAR(rc.y,2,1e-12); EXPECT_NEAR(rc.z,2,1e-12);
    cout << "       cubic = " << cubic << endl;

    UnitCell ortho(4,5,6,90,90,90);
    EXPECT_NEAR(ortho.GetCellVolume(),120,1e-9);
    EXPECT_NEAR(ortho.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(ortho.GetDistance(rvec3_t{1,1,1}),sqrt(4*4+5*5+6*6),1e-9);
    EXPECT_EQ(ortho.GetNumCells(9),Vector3D<int>(3,2,2));
    rvec3_t ro=ortho.ToCartesian({1,1,1});
    EXPECT_NEAR(ro.x,4,1e-12); EXPECT_NEAR(ro.y,5,1e-12); EXPECT_NEAR(ro.z,6,1e-12);
    cout << "orthorhombic = " << ortho << endl;

    UnitCell tri(4,5,6,80,95,75);
    EXPECT_NEAR(tri.GetCellVolume(),113.04412274615466,1e-9);
    EXPECT_NEAR(tri.GetMinimumCellEdge(),4,1e-9);
    EXPECT_NEAR(tri.GetDistance(rvec3_t{1,1,1}),9.6740982428456377,1e-9);
    EXPECT_EQ(tri.GetNumCells(12),Vector3D<int>(4,3,3)); //interplanar spacing, not edge length
    cout << "   triclinic = " << tri << endl;

    //  Reciprocal cell B = 2π A⁻ᵀ: volume is (2π)³/V, and reciprocal-of-
    //  reciprocal returns the original cell (the 2π convention round-trips).
    UnitCell recip=tri.MakeReciprocalCell();
    EXPECT_NEAR(recip.GetCellVolume(),pow(2*Pi,3)/tri.GetCellVolume(),1e-12);
    UnitCell back=recip.MakeReciprocalCell();
    EXPECT_NEAR(back.GetCellVolume(),tri.GetCellVolume(),1e-9);
    EXPECT_NEAR(back.GetDistance({1,0,0}),tri.GetDistance({1,0,0}),1e-9);
    EXPECT_NEAR(back.GetDistance({1,1,1}),tri.GetDistance({1,1,1}),1e-9);
}

TEST_F(StructureTests, Lattice)
{
    double a_0=0.529177; //Ångstrom
    UnitCell SiCell(5.43/a_0); //Convert Ångstrom to atomic units a.u.
    SiCell.Insert(new Atom(14,0,rvec3_t{ 0, 0, 0}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.5,.5, 0}));
    SiCell.Insert(new Atom(14,0,rvec3_t{ 0,.5,.5}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.5, 0,.5}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.25,.25,.25}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.75,.75,.25}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.25,.75,.75}));
    SiCell.Insert(new Atom(14,0,rvec3_t{.75,.25,.75}));
    Lattice Si(SiCell,ivec3_t(2,2,2));
    // cout << "Si lattice = " << Si << endl;
    EXPECT_EQ(SiCell.GetNumAtoms(),8);
    EXPECT_EQ(SiCell.GetNuclearCharge(),8*14);
    EXPECT_EQ(SiCell.GetNetCharge(),0);
    EXPECT_EQ(SiCell.GetNumElectrons(),8*14);
    EXPECT_NEAR(Si.GetLatticeVolume(),pow(2*5.43/a_0,3),1e-6); //2x2x2 supercell = (2*edge)^3
    EXPECT_EQ(Si.GetNumSites(),64);
    EXPECT_EQ(Si.GetNumBasisSites(),8);
    EXPECT_EQ(Si.GetNumUnitCells(),8);
    EXPECT_EQ(Si.GetSiteNumber (rvec3_t(1.75,0.75,1.25)),45);
    EXPECT_EQ(Si.GetBasisNumber(rvec3_t(1.75,0.75,1.25)),5);
    EXPECT_EQ(Si.GetBasisNumber(13),5);

    for (size_t i=0;i<Si.GetNumSites();i++)
        EXPECT_EQ(Si.GetSiteNumber(Si.GetCoordinate(i)),i);
 
    double Emax=2.0;
    Lattice Rl=Si.Reciprocal(Emax);
    // cout << "Si reciprocal lattice = " << Rl << endl;
    std::vector<ivec3_t> cells =Rl.GetCellsInSphere(Emax);
    UnitCell Rc=Rl.GetUnitCell();
    for (auto c:cells)
        EXPECT_LT(Rc.GetDistance(c),Emax);

    std::vector<rvec3_t> ks=Si.GetReciprocalGrid();
    cout << "BZ grid:" << endl;
    for (auto k:ks) cout << "   "  << k << endl;
}