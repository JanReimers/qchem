// File: UnitTests/BasisSet_Lattice.C  Unit test the Lattice IBS Evaluators and basis sets.

#include "gtest/gtest.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

class BasisSet_Lattice : public ::testing::Test
{
    BasisSet_Lattice()
    : N(2)
    , SiBasis(new Molecule())
    , SiCell(5.43/a_0) //Convert lattice constant to atom units (au).
    , Si(SiCell,ivec3_t(N,N,N),std::shared_ptr<Cluster>(SiBasis))
    {
        SiBasis->Insert(new Atom(14,0,rvec3_t{ 0, 0, 0}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{.5,.5, 0}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{ 0,.5,.5}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{.5, 0,.5}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{.25,.25,.25}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{.75,.75,.25}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{.25,.75,.75}));
        SiBasis->Insert(new Atom(14,0,rvec3_t{.75,.25,.75}));
    }

    const double a_0=0.529177; //Ångstrom per A.U.
    size_t    N;
    Molecule* SiBasis;    
    UnitCell  SiCell;
    Lattice   Si;
};

TESTF(BasisSet_Lattice,ConstructBasisSet)
{
    BasisSet bs=new PlaneWave::BasisSet(Si,10.0);
    cout << *bs << endl;
}