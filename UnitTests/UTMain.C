// UTMain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    

int main(int argc, char **argv)
{
     testing::InitGoogleTest(&argc, argv);
//     testing::GTEST_FLAG(filter) = "A_PG_HF_P_92*";

//     testing::GTEST_FLAG(filter) = "MeshIntegralsTests.GObritals:M_PG_*:Multiple/A_PG*";
//  testing::GTEST_FLAG(filter) = "Multiple/A_PG_DFT*";
    testing::GTEST_FLAG(filter) = "Multiple/A_SLm_HF_P.Multiple/2";
//
//    testing::GTEST_FLAG(filter) = "STLTesting.RangeBasedLoops";
//      testing::GTEST_FLAG(filter) = "GaussianRadialERITests.*";
//      testing::GTEST_FLAG(filter) = "GaussianRadialIntegralTests.*";
//      testing::GTEST_FLAG(filter) = "SlaterRadialIntegralTests.*";
//    testing::GTEST_FLAG(filter) = "PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_EigenSolverTests.*";

    return RUN_ALL_TESTS();
}

