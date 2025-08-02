#include "gtest/gtest.h"    

using std::cout;
using std::endl;

import qchem.BasisSet.Lattice.PlaneWave;
// import qchem.Atom;
// import Cluster.UnitCell;
// import qchem.Lattice;
// import qchem.Molecule;

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

Molecule* MakeSiBasis()
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
    return SiBasis;
}

class BasisSetTests : public ::testing::Test
{
    BasisSetTests() 
        : SiCell(5.43/a_0) //Convert Ångstrom to atomic units a.u.
        , Si(SiCell,IVec3(2,2,2),std::shared_ptr<Cluster>(MakeSiBasis()))
        {}

    static const double a_0=0.529177; //Ångstrom
    
    UnitCell SiCell; 
    Lattice  Si;
}; 

TEST_F(BasisSetTests, PlaneWave)
{
    // BasisSet* bs=new PlaneWave::BasisSet(Si,2.0);
}