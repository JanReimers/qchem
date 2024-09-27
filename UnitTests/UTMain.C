// UTMain.C Main entry point for all unit tests.

#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

//    testing::GTEST_FLAG(filter) = "SlaterIntegralTests.*";
//    testing::GTEST_FLAG(filter) = "AtomsHFPolarized/HartreeFockAtomTester.AtomsHFPolarized/1";
//    testing::GTEST_FLAG(filter) = "SemiHartreeFockMoleculeTester.*";
//    testing::GTEST_FLAG(filter) = "DFTMoleculeTester.N2";
    testing::GTEST_FLAG(filter) = "STLTesting.RangeBasedLoops";

    
    return RUN_ALL_TESTS();
}

