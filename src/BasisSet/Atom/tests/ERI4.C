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

TEST_F(ERI4Tests,MatMulTimings)
{
#ifdef DEBUG
    GTEST_SKIP() << "Timing is unreliable in Debug; runs in Release only.";
#endif
    size_t Nrep=100;
    size_t Nab=60,Ncd=60;
    ERI4 Jabcd(Nab,Ncd);
    rsmat_t Dab(Nab),Dcd(Ncd);
    random(Jabcd);
    random(Dab);
    random(Dcd);
    std::chrono::duration<double> elapsed_seconds1,elapsed_seconds2;
    {
        rsmat_t Jab=blazem::zero<double>(Nab);
        auto start = std::chrono::system_clock::now();
        for (size_t i=0;i<Nrep;i++)
            MatMul(Jab,Jabcd,Dcd);
        auto end = std::chrono::system_clock::now();
        elapsed_seconds1 = end-start;
        std::cout << "MatMul(Jabcd,Dcd) elapsed time: " << elapsed_seconds1.count() << "s" << std::endl;
    }
    {        
        rsmat_t Jcd=blazem::zero<double>(Ncd);
        auto start = std::chrono::system_clock::now();
        for (size_t i=0;i<Nrep;i++)
            MatMul(Jcd,Dab,Jabcd);
        auto end = std::chrono::system_clock::now();
        elapsed_seconds2 = end-start;
        std::cout << "MatMul(Dab,Jabcd) elapsed time: " << elapsed_seconds2.count() << "s" << std::endl;
    }
    EXPECT_GT(elapsed_seconds2,elapsed_seconds1);
#ifndef DEBUG
    //Ratio is much lower in debug mode.
    //Ration is also much when running the full testsuite.
    EXPECT_GT(elapsed_seconds2/elapsed_seconds1,2.5); 
#endif
}

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