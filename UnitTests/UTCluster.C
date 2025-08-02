#include <cmath>
#include <memory>
#include "gtest/gtest.h"    

using std::cout;
using std::endl;

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
    cout << "default atom=" << atom << endl;
    cout << "           B=" << B << endl;
    cout << "          F-=" << F_minus << endl;
    cout << "         Na+=" << Na_plus << endl;
  
}

TEST_F(ClusterTests, UnitCell)
{
    UnitCell cubic(4,4,4,90,90,90);
    EXPECT_EQ(cubic.GetCellVolume(),64);
    EXPECT_EQ(cubic.GetMinimumCellEdge(),4);
    EXPECT_EQ(cubic.GetDistance(RVec3{1,1,1}),4*sqrt(3));
    EXPECT_EQ(cubic.GetNumCells(9),Vector3D<int>(3,3,3));
    cout << "       cubic = " << cubic << endl;
    UnitCell ortho(4,5,6,90,90,90);
    EXPECT_EQ(ortho.GetCellVolume(),120);
    EXPECT_EQ(ortho.GetMinimumCellEdge(),4);
    EXPECT_EQ(ortho.GetDistance(RVec3{1,1,1}),sqrt(4*4+5*5+6*6));
    EXPECT_EQ(ortho.GetNumCells(9),Vector3D<int>(3,2,2));
    cout << "orthorhombic = " << ortho << endl;
    UnitCell tri(4,5,6,80,95,75);
    EXPECT_EQ(tri.GetCellVolume(),113.04412274615466);
    EXPECT_EQ(tri.GetMinimumCellEdge(),4);
    EXPECT_EQ(tri.GetDistance(RVec3{1,1,1}),9.6740982428456377);
    EXPECT_EQ(tri.GetNumCells(12),Vector3D<int>(3,3,2));
    cout << "   triclinic = " << tri << endl;
}

TEST_F(ClusterTests, Lattice)
{
    double a_0=0.529177; //Ångstrom
    Molecule* SiBasis=new Molecule();
    SiBasis->Insert(new Atom(14,0,RVec3{ 0, 0, 0}));
    SiBasis->Insert(new Atom(14,0,RVec3{.5,.5, 0}));
    SiBasis->Insert(new Atom(14,0,RVec3{ 0,.5,.5}));
    SiBasis->Insert(new Atom(14,0,RVec3{.5, 0,.5}));
    SiBasis->Insert(new Atom(14,0,RVec3{.25,.25,.25}));
    SiBasis->Insert(new Atom(14,0,RVec3{.75,.75,.25}));
    SiBasis->Insert(new Atom(14,0,RVec3{.25,.75,.75}));
    SiBasis->Insert(new Atom(14,0,RVec3{.75,.25,.75}));
    UnitCell SiCell(5.43/a_0); //Convert Ångstrom to atomic units a.u.
    Lattice Si(SiCell,IVec3(2,2,2),std::shared_ptr<Cluster>(SiBasis));
    // cout << "Si lattice = " << Si << endl;
    EXPECT_EQ(Si.GetNumAtoms(),8);
    EXPECT_EQ(Si.GetNuclearCharge(),8*14);
    EXPECT_EQ(Si.GetNetCharge(),0);
    EXPECT_EQ(Si.GetNumElectrons(),8*14);
    EXPECT_EQ(Si.GetLatticeVolume(),pow(4*5.43/a_0,3));
    EXPECT_EQ(Si.GetNumSites(),64);
    EXPECT_EQ(Si.GetNumBasisSites(),8);
    EXPECT_EQ(Si.GetNumUnitCells(),8);
    EXPECT_EQ(Si.GetSiteNumber (RVec3(1.75,0.75,1.25)),45);
    EXPECT_EQ(Si.GetBasisNumber(RVec3(1.75,0.75,1.25)),5);
    EXPECT_EQ(Si.GetBasisNumber(13),5);

    for (size_t i=0;i<Si.GetNumSites();i++)
        EXPECT_EQ(Si.GetSiteNumber(Si.GetCoordinate(i)),i);
 
    double Emax=2.0;
    Lattice Rl=Si.Reciprocal(Emax);
    // cout << "Si reciprocal lattice = " << Rl << endl;
    std::vector<IVec3> cells =Rl.GetCellsInSphere(Emax);
    UnitCell Rc=Rl.GetUnitCell();
    for (auto c:cells)
        EXPECT_LT(Rc.GetDistance(c),Emax);

    std::vector<RVec3> ks=Si.GetReciprocalGrid();
    cout << "BZ grid:" << endl;
    for (auto k:ks) cout << "   "  << k << endl;
}