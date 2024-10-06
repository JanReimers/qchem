// UTMain.C Main entry point for all unit tests.

#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

//    testing::GTEST_FLAG(filter) = "SlaterIntegralTests.*";
//    testing::GTEST_FLAG(filter) = "PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_EigenSolverTests.*";

//    testing::GTEST_FLAG(filter) = "SemiHartreeFockMoleculeTester.*";
//    testing::GTEST_FLAG(filter) = "DFTMoleculeTester.*"; 
//    testing::GTEST_FLAG(filter) = "DFTMoleculeTester.*:MoleculesDFTPolarized/DFTMoleculeTester.*";
//    testing::GTEST_FLAG(filter) = "STLTesting.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "HartreeFockMoleculeTester.*";

//    testing::GTEST_FLAG(filter) = "HartreeFockAtomTester.AtomsHFEigenSolvers";
//    testing::GTEST_FLAG(filter) = "AtomsSemiDFTPolarized/SemiHartreeFockAtomTester.AtomsSemiDFTPolarized/*";
//    testing::GTEST_FLAG(filter) = "AtomsDFTPolarized/DFTAtomTester.*";
//    testing::GTEST_FLAG(filter) = "AtomsHFPolarized/HartreeFockAtomTester.AtomsHFPolarized/*";
    return RUN_ALL_TESTS();
}

