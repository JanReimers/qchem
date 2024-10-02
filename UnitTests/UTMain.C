// UTMain.C Main entry point for all unit tests.

#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

//    testing::GTEST_FLAG(filter) = "SlaterIntegralTests.*";
//        testing::GTEST_FLAG(filter) = "AtomsHFPolarized/HartreeFockAtomTester.AtomsHFPolarized/1";
//    testing::GTEST_FLAG(filter) = "SemiHartreeFockMoleculeTester.*";
//    testing::GTEST_FLAG(filter) = "DFTMoleculeTester.N2";
//    testing::GTEST_FLAG(filter) = "STLTesting.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "AtomsHFPolarized/HartreeFockAtomTester.AtomsHFPolarized/2";
//    testing::GTEST_FLAG(filter) = "MoleculesHFPolarized/HartreeFockMoleculeTester.MoleculesHFPolarized/*";
//    testing::GTEST_FLAG(filter) = "AtomsDFTPolarized/DFTAtomTester.AtomsDFTPolarized/0";
//    testing::GTEST_FLAG(filter) = "PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_EigenSolverTests.*";
//    testing::GTEST_FLAG(filter) = "HartreeFockAtomTester.AtomsHFEigenSolvers";
    
    return RUN_ALL_TESTS();
}

