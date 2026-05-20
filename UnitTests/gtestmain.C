// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    
import qchem.BasisSet.Internal.DB_Cache_RAM;

int main(int argc, char **argv)
{
     BasisSet::theGlobalCache=new BasisSet::IntegralsCache_RAM<double>(true);  
     testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}

