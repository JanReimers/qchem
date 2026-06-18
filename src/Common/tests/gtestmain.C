// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    
// import qchem.BasisSet.Internal.DB_Cache_RAM;


// class GlobalEnvironment : public ::testing::Environment {
//  public:
//   ~GlobalEnvironment() override {}

//   // Called once before all tests run
//   void SetUp() override {
//     std::cout << "Global Environment Setup" << std::endl;
//     BasisSet::theGlobalCache=new BasisSet::IntegralsCache_RAM<double>(true);  
//   }

//   // Called once after all tests finish
//   void TearDown() override {
//     std::cout << "Global Environment TearDown" << std::endl;
//     delete BasisSet::theGlobalCache;  
//   }
// };


int main(int argc, char **argv)
{
     
     testing::InitGoogleTest(&argc, argv);
     // Register the global environment
    //  ::testing::AddGlobalTestEnvironment(new GlobalEnvironment);
     return RUN_ALL_TESTS();
}


