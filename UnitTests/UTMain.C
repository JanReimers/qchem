// UTMain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    

int main(int argc, char **argv)
{
     testing::InitGoogleTest(&argc, argv);
     // testing::GTEST_FLAG(filter) = "SlaterRadialIntegralTests.*";  
     // testing::GTEST_FLAG(filter) = "GaussianRadialIntegralTests.*";
     // testing::GTEST_FLAG(filter) = "DiracIntegralTests.*";  
//     testing::GTEST_FLAG(filter) = "DiracIntegralTests.*:Multiple/A_*_DHF.Multiple/*";
     //    testing::GTEST_FLAG(filter) = "Multiple/A_*_DHF.Multiple/*:DiracIntegralTests.*";
     testing::GTEST_FLAG(filter) = "Multiple/A_SL*_*_*.Multiple/*:SlaterRadialIntegralTests.*";
//     testing::GTEST_FLAG(filter) = "Multiple/A_SLmj_DHF.Multiple/*";
     // testing::GTEST_FLAG(filter) = "Multiple/A_*_DHF.Multiple/*";
//     testing::GTEST_FLAG(filter) = "A_SG_DHF_H.Phir";
//     testing::GTEST_FLAG(filter) = "A_SG_HFP_H.Phir";
     //   testing::GTEST_FLAG(filter) = "Multiple/A_*_DFT_P.Multiple/*:Multiple/A_SLmj_HF.Multiple/*";

//     testing::GTEST_FLAG(filter) = "libCintTests.*"; 
//     testing::GTEST_FLAG(filter) = "ERI4Tests.*"; 
//     testing::GTEST_FLAG(filter) = "Multiple/A_SG_DFT_P.Multiple/*";
//     testing::GTEST_FLAG(filter) = "Multiple/A_SGm_HF_P.Multiple/*:Multiple/A_SG_HF_P.Multiple/*:Multiple/A_SLm_HF_P.Multiple/*:Multiple/A_SL_HF_P.Multiple/*";
//     testing::GTEST_FLAG(filter) = "Multiple/A_SGm_HF_P.Multiple/*:Multiple/A_SG_HF_P.Multiple/*";
//     testing::GTEST_FLAG(filter) = "Multiple/A_SGm_HF_P.Multiple/*";
     //   testing::GTEST_FLAG(filter) = "Multiple/A_PG_HF*";
     // testing::GTEST_FLAG(filter) = "M_*";

//     testing::GTEST_FLAG(filter) = "MeshIntegralsTests.GObritals:M_PG_*:Multiple/A_PG*";
//      testing::GTEST_FLAG(filter) = "Multiple/A_SG_HF_*.Multiple/*";
//    testing::GTEST_FLAG(filter) = "STLTesting.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "GaussianRadialERITests.*";
//      testing::GTEST_FLAG(filter) = "SlaterRadialIntegralTests.*";
//    testing::GTEST_FLAG(filter) = "PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "qchem_PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "ElectronConfigurationTests.*";

    return RUN_ALL_TESTS();
}

