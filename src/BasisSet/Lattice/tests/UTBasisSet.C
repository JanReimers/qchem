#include "gtest/gtest.h"    

using std::cout;
using std::endl;

// import qchem.BasisSet.Lattice.PlaneWave;
import qchem.Structure;
import qchem.UnitCell;
import qchem.Lattice;


class BasisSetTests : public ::testing::Test
{
public: 
    BasisSetTests() 
        : SiCell(5.43/a_0) //Convert Ångstrom to atomic units a.u.
        , Si(SiCell,ivec3_t(2,2,2))
        {
            SiCell.Insert(new Atom(14,0,rvec3_t{ 0, 0, 0}));
            SiCell.Insert(new Atom(14,0,rvec3_t{.5,.5, 0}));
            SiCell.Insert(new Atom(14,0,rvec3_t{ 0,.5,.5}));
            SiCell.Insert(new Atom(14,0,rvec3_t{.5, 0,.5}));
            SiCell.Insert(new Atom(14,0,rvec3_t{.25,.25,.25}));
            SiCell.Insert(new Atom(14,0,rvec3_t{.75,.75,.25}));
            SiCell.Insert(new Atom(14,0,rvec3_t{.25,.75,.75}));
            SiCell.Insert(new Atom(14,0,rvec3_t{.75,.25,.75}));
        }

    constexpr static const double a_0=0.529177; //Ångstrom
    
    UnitCell SiCell; 
    Lattice  Si;
}; 

TEST_F(BasisSetTests, PlaneWave)
{
    // BasisSet* bs=new PlaneWave::BasisSet(Si,2.0);
}