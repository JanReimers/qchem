// File: UnitTests/ERI4.C  Test the four index supermatrix used for storing 2 electron repulsion (ERI) integrals

#include "gtest/gtest.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include <blaze/math/expressions/DMatDMatEqualExpr.h> //op== inside gtest header.
using std::cout;
using std::endl;
import qchem.BasisSet.Internal.ERI4;
import qchem.Blaze;
using namespace qchem;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class ERI4Tests : public ::testing::Test
{

};

void MatMul(rsmat_t& Scd, const rsmat_t& Sab, const ERI4& gabcd)
{
    size_t Nab=gabcd.Nab();
    assert(Scd.rows()==gabcd(0,0).rows());
    for (auto ia:iv_t(0,Nab))
    {
        Scd+=gabcd(ia,ia)*Sab(ia,ia);
        for (auto ib:iv_t(ia+1,Nab))
            Scd+=2*gabcd(ia,ib)*Sab(ia,ib);
    }
}


void random(rsmat_t& s)
{
    for (auto i:iv_t(0,s.rows()))
        for (auto j:iv_t(i,s.rows()))
            s(i,j)=random();
}

void random(ERI4& Jabcd)
{
    for (auto a:iv_t(0,Jabcd.Nab()))
        for (auto b:iv_t(a,Jabcd.Nab()))
            random(Jabcd(a,b));
}

// MatMulTimings DELETED 2026-08-12 (user).  It asserted a wall-clock RATIO
// (elapsed(Dab*Jabcd) / elapsed(Jabcd*Dcd) > 2.5), which is not a property of the code but of the machine
// it runs on: it failed spuriously whenever the box was busy -- its own comment already conceded "the ratio
// is also much [lower] when running the full testsuite" -- and it skipped itself entirely in Debug.  A test
// that has to be switched off in one build and is unreliable in the other is measuring the wrong thing.
//
// The FACT it recorded is worth keeping, so it moved to the ERI4::MatMul declaration where a reader will
// meet it: contracting over the cd index is several times faster than contracting over ab, because the ab
// side is the outer (block) index and the cd side is contiguous within a block.  That belongs in a
// benchmark (src/BasisSet/Molecule/bench) if it is ever to be MEASURED again, not in a correctness suite.

TEST_F(ERI4Tests,Transpose)
{
    size_t Nab=30,Ncd=40;
    ERI4 Jabcd(Nab,Ncd);
    random(Jabcd);
    ERI4 Jcdab=Jabcd.Transpose();
    rsmat_t Dab(Nab),Dcd(Ncd);
    random(Dab);
    random(Dcd);

    rsmat_t Jab1=blazem::zero<double>(Nab), Jab2=blazem::zero<double>(Nab);
    MatMul(Jab1,Jabcd,Dcd);
    MatMul(Jab2,Dcd,Jcdab);
    EXPECT_EQ(Jab1,Jab2);
    
    rsmat_t Jcd1=blazem::zero<double>(Ncd), Jcd2=blazem::zero<double>(Ncd);
    MatMul(Jcd1,Dab,Jabcd);
    MatMul(Jcd2,Jcdab,Dab);
    EXPECT_EQ(Jcd1,Jcd2);
}

// The fused ScatterBoth must reproduce the TWO independent contractions it replaces in the Fock build:
//   Si += J·Dj   (localized, == free MatMul(Si,J,Dj))
//   Sj += Jᵀ·Di  (== the whole-block-add MatMul(Sj,Di,J) above -- same order & weights, so bit-identical).
TEST_F(ERI4Tests,ScatterBoth)
{
    size_t Nab=30,Ncd=40;
    ERI4 Jabcd(Nab,Ncd);
    random(Jabcd);
    rsmat_t Di(Nab),Dj(Ncd);
    random(Di);
    random(Dj);

    rsmat_t Si_ref=blazem::zero<double>(Nab), Sj_ref=blazem::zero<double>(Ncd);
    MatMul(Si_ref,Jabcd,Dj);   // Si += J·Dj  (production free MatMul)
    MatMul(Sj_ref,Di,Jabcd);   // Sj += Jᵀ·Di (test-local whole-block-add MatMul)

    rsmat_t Si=blazem::zero<double>(Nab), Sj=blazem::zero<double>(Ncd);
    Jabcd.ScatterBoth(Si,Sj,Di,Dj);   // one pass, both targets

    EXPECT_EQ(Si,Si_ref);
    EXPECT_EQ(Sj,Sj_ref);
}