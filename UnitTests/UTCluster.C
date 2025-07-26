#include <cmath>
#include "gtest/gtest.h"    

import qchem.Atom;
import Cluster.UnitCell;
import qchem.Lattice;
import qchem.Molecule;

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

class ClusterTests : public ::testing::Test
{

}; 

TEST_F(ClusterTests, Atom)
{
    Atom atom;    
    Atom B(5,0);
    Atom F_minus(9,-1);
    Atom Na_plus(11,1,{1,1,1});
    EXPECT_EQ(atom.GetNumElectrons(),0);
    EXPECT_EQ(B.GetNumElectrons(),5);
    EXPECT_EQ(F_minus.GetNumElectrons(),10);
    EXPECT_EQ(Na_plus.GetNumElectrons(),10);
    std::cout << "default atom=" << atom << std::endl;
    std::cout << "           B=" << B << std::endl;
    std::cout << "          F-=" << F_minus << std::endl;
    std::cout << "         Na+=" << Na_plus << std::endl;
  
}

TEST_F(ClusterTests, UnitCell)
{
    UnitCell cubic(4,4,4,90,90,90);
    EXPECT_EQ(cubic.GetCellVolume(),64);
    EXPECT_EQ(cubic.GetMinimumCellEdge(),4);
    EXPECT_EQ(cubic.GetDistance(RVec3{1,1,1}),4*sqrt(3));
    EXPECT_EQ(cubic.GetNumCells(9),Vector3D<int>(3,3,3));
    std::cout << "       cubic = " << cubic << std::endl;
    UnitCell ortho(4,5,6,90,90,90);
    EXPECT_EQ(ortho.GetCellVolume(),120);
    EXPECT_EQ(ortho.GetMinimumCellEdge(),4);
    EXPECT_EQ(ortho.GetDistance(RVec3{1,1,1}),sqrt(4*4+5*5+6*6));
    EXPECT_EQ(ortho.GetNumCells(9),Vector3D<int>(3,2,2));
    std::cout << "orthorhombic = " << ortho << std::endl;
    UnitCell tri(4,5,6,80,95,75);
    EXPECT_EQ(tri.GetCellVolume(),113.04412274615466);
    EXPECT_EQ(tri.GetMinimumCellEdge(),4);
    EXPECT_EQ(tri.GetDistance(RVec3{1,1,1}),9.6740982428456377);
    EXPECT_EQ(tri.GetNumCells(12),Vector3D<int>(3,3,2));
    std::cout << "   triclinic = " << tri << std::endl;
}

TEST_F(ClusterTests, Lattice)
{
    Molecule* SiBasis=new Molecule();
    SiBasis->Insert(new Atom(14,0,RVec3{ 0, 0, 0}));
    SiBasis->Insert(new Atom(14,0,RVec3{.5,.5, 0}));
    SiBasis->Insert(new Atom(14,0,RVec3{ 0,.5,.5}));
    SiBasis->Insert(new Atom(14,0,RVec3{.5, 0,.5}));
    SiBasis->Insert(new Atom(14,0,RVec3{.25,.25,.25}));
    SiBasis->Insert(new Atom(14,0,RVec3{.75,.75,.25}));
    SiBasis->Insert(new Atom(14,0,RVec3{.25,.75,.75}));
    SiBasis->Insert(new Atom(14,0,RVec3{.75,.25,.75}));
    UnitCell SiCell(5.43);
    Lattice Si(SiCell,Vector3D<int>(1,1,1),SiBasis);
    std::cout << "Si lattice = " << Si << std::endl;
}