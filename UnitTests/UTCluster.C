#include <cmath>
#include "gtest/gtest.h"    

import qchem.Atom;
import Cluster.UnitCell;

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
    std::cout << "cubic = " << cubic << std::endl;
}