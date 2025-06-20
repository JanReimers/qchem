// UTMain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    

int main(int argc, char **argv)
{
     testing::InitGoogleTest(&argc, argv);
     //testing::GTEST_FLAG(filter) = "A_BS_1E_U.*";  
     // testing::GTEST_FLAG(filter) = "Multiple/A_SL_HF_P.Multiple/*";  
     // testing::GTEST_FLAG(filter) = "Multiple/A_BS_HF_U.Multiple/*";  
     // testing::GTEST_FLAG(filter) = "Multiple/A_BS*_HF_*.Multiple/*";  
     // testing::GTEST_FLAG(filter) = "BSplineTests.*";   
     // testing::GTEST_FLAG(filter) = "SymQNTests.*";  
     // testing::GTEST_FLAG(filter) = "SlaterRadialIntegralTests.*";  
     // testing::GTEST_FLAG(filter) = "GaussianRadialIntegralTests.*";
     // testing::GTEST_FLAG(filter) = "DiracIntegralTests.*";  
//     testing::GTEST_FLAG(filter) = "DiracIntegralTests.*:Multiple/A_*_DHF.Multiple/*";
     //    testing::GTEST_FLAG(filter) = "Multiple/A_*_DHF.Multiple/*:DiracIntegralTests.*";
     // testing::GTEST_FLAG(filter) = "M_PG_*_*.*";
     //  testing::GTEST_FLAG(filter) = "ElectronConfigurationTests.*";
     //  testing::GTEST_FLAG(filter) = "OrthogonalizeTests.*";
     //   testing::GTEST_FLAG(filter) = "jsonTests.*";
//     testing::GTEST_FLAG(filter) = "Multiple/A_SLmj_DHF.Multiple/*";
     // testing::GTEST_FLAG(filter) = "Multiple/A_PG_HF_P.Multiple/1";
//     testing::GTEST_FLAG(filter) = "A_SG_DHF_H.Phir";
     //  testing::GTEST_FLAG(filter) = "A_SG_HFP_H.Phir";
     //   testing::GTEST_FLAG(filter) = "Multiple/*_*_DFT_*.Multiple/*";
     //   testing::GTEST_FLAG(filter) = "Multiple/A_SG_DFT_U.Multiple/2";

//     testing::GTEST_FLAG(filter) = "libCintTests.*"; 
//     testing::GTEST_FLAG(filter) = "ERI4Tests.*"; 
//    testing::GTEST_FLAG(filter) = "STLTesting.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "GaussianRadialERITests.*";
//      testing::GTEST_FLAG(filter) = "SlaterRadialIntegralTests.*";
//    testing::GTEST_FLAG(filter) = "PersistanceTests.*";
//    testing::GTEST_FLAG(filter) = "OrthogonalizeTests.*";

    return RUN_ALL_TESTS();
}

