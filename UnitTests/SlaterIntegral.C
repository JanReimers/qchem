// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Wigner3j.H"
#include "Misc/DFTDefines.H"
#include "oml/imp/ran250.h"
#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class SlaterIntegralTests : public ::testing::Test
{
public:
    SlaterIntegralTests() :ab(1.5), cd(0.1) {};
    void Randomize()
    {
        ab=1000.*exp(-OMLRandPos<double>()*10.0);
        cd=1000.*exp(-OMLRandPos<double>()*10.0);
    }

    double ab,cd;
};

TEST_F(SlaterIntegralTests, R0_0000)
{
    SlaterIntegrals R(ab,cd);
    double expected=Pi12/8*(1/(ab*cd*sqrt(ab+cd)));
    EXPECT_DOUBLE_EQ(expected,R(0,0,0,0,0));
}

TEST_F(SlaterIntegralTests, R0_0000_Random)
{
    for (int i=0; i<1000; i++)
    {
        Randomize();
        SlaterIntegrals R(ab,cd);
        double expected=Pi12/8*(1/(ab*cd*sqrt(ab+cd)));
        EXPECT_DOUBLE_EQ(expected,R(0,0,0,0,0));
    }
}

TEST_F(SlaterIntegralTests, R0_1100)
{
    SlaterIntegrals R(ab,cd);
    double expected=1/ab+1/(2*(ab+cd));
    EXPECT_DOUBLE_EQ(expected,R(0,1,1,0,0)/R(0,0,0,0,0));
}

TEST_F(SlaterIntegralTests, R0_1100_Random)
{
    for (int i=0; i<1000; i++)
    {
        Randomize();
        SlaterIntegrals R(ab,cd);
        double expected=R(0,0,0,0,0)*(1/ab+1/(2*(ab+cd)));
        EXPECT_DOUBLE_EQ(expected,R(0,1,1,0,0));
    }
}

TEST_F(SlaterIntegralTests, R0_0011)
{
    SlaterIntegrals R(ab,cd);
    double expected=R(0,0,0,0,0)*(1/cd+1/(2*(ab+cd)));
    EXPECT_DOUBLE_EQ(expected,R(0,0,0,1,1));
}

TEST_F(SlaterIntegralTests, R0_0011_Random)
{
    for (int i=0; i<1000; i++)
    {
        Randomize();
        SlaterIntegrals R(ab,cd);
        double expected=R(0,0,0,0,0)*(1/cd+1/(2*(ab+cd)));
        EXPECT_DOUBLE_EQ(expected,R(0,0,0,1,1));
    }
}
//
TEST_F(SlaterIntegralTests, R1_1010)
{
    SlaterIntegrals R(ab,cd);
    double expected=R(0,0,0,0,0)*(3/(2*(ab+cd)));
    EXPECT_DOUBLE_EQ(expected,R(1,1,0,1,0));
}

TEST_F(SlaterIntegralTests, R0_1111)
{
    SlaterIntegrals R(ab,cd);
    double expected=R(0,0,0,0,0)*3/4*(2/(ab*cd)+1/((ab+cd)*(ab+cd)));
    EXPECT_DOUBLE_EQ(expected,R(0,1,1,1,1));
}
TEST_F(SlaterIntegralTests, R2_1111)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd2=abcd*abcd;
    double expected=R(0,0,0,0,0)*15/4*(1/abcd2);
    EXPECT_DOUBLE_EQ(expected,R(2,1,1,1,1));
}

TEST_F(SlaterIntegralTests, R2_0202)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd2=abcd*abcd;
    double expected=R(0,0,0,0,0)*15/4*(1/abcd2);
    EXPECT_DOUBLE_EQ(expected,R(2,0,2,0,2));
}

TEST_F(SlaterIntegralTests, R1_1212)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd3=abcd*abcd*abcd;
    double expected=15/(4*abcd3)*(cd/ab+ab/cd+7./2.);
    EXPECT_DOUBLE_EQ(expected,R(1,1,2,1,2)/R(0,0,0,0,0));
}

TEST_F(SlaterIntegralTests, R3_1212)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd3=abcd*abcd*abcd;
    double expected=R(0,0,0,0,0)*105/(8*abcd3);
    EXPECT_DOUBLE_EQ(expected,R(3,1,2,1,2));
}
//
TEST_F(SlaterIntegralTests, R4_2222)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd4=abcd*abcd*abcd*abcd;
    double expected=R(0,0,0,0,0)*945/(16*abcd4);
    EXPECT_DOUBLE_EQ(expected,R(4,2,2,2,2));
}
//
TEST_F(SlaterIntegralTests, R2_2222)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd4=abcd*abcd*abcd*abcd;
    double expected=R(0,0,0,0,0)*105/(8*abcd4)*(cd/ab+ab/cd+9./2.);
    EXPECT_DOUBLE_EQ(expected,R(2,2,2,2,2));
}
//
//// This one will be fun!!!
TEST_F(SlaterIntegralTests, R0_2222)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd4=abcd*abcd*abcd*abcd;
    double expected=15/(4*abcd4)*(2*abcd*(cd/(ab*ab)+ab/(cd*cd)) + 7*(cd/ab+ab/cd)+63./4.);
    EXPECT_DOUBLE_EQ(expected,R(0,2,2,2,2)/R(0,0,0,0,0));
}

//// This one will be more fun!!!
TEST_F(SlaterIntegralTests, R0_1122)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double abcd3=abcd*abcd*abcd;
    double Iab=7/(2*abcd)+1/ab;
    double Icd=5/(2*abcd)+1/cd;
    double Icd2=5/(2*abcd*abcd)+1/(cd*cd);
    double expected=3/(2*abcd3)*(5*cd/2*Iab+ab*abcd*(Icd*Icd+Icd2));
    EXPECT_DOUBLE_EQ(expected,R(0,1,1,2,2)/R(0,0,0,0,0));
}

//// This one should be easier
TEST_F(SlaterIntegralTests, R0_0022)
{
    SlaterIntegrals R(ab,cd);
    double abcd=ab+cd;
    double expected=1/(abcd)*((2*ab+3*cd)/(cd*cd)+3/(4*abcd));
    EXPECT_DOUBLE_EQ(expected,R(0,0,0,2,2)/R(0,0,0,0,0));
}

TEST_F(SlaterIntegralTests, Wigner1)
{
    Wigner3j w;
    double expected=3./70.;
    EXPECT_DOUBLE_EQ(w(1,2,3),expected);
    EXPECT_DOUBLE_EQ(w(1,3,2),expected);
    EXPECT_DOUBLE_EQ(w(2,1,3),expected);
    EXPECT_DOUBLE_EQ(w(2,3,1),expected);
    EXPECT_DOUBLE_EQ(w(3,2,1),expected);
    EXPECT_DOUBLE_EQ(w(3,1,2),expected);
}
