// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    
import qchem.BasisSet.DB_Cache;

int main(int argc, char **argv)
{
      if (BasisSet::theGlobalCache==0)
          BasisSet::theGlobalCache=new BasisSet::IntegralsCache_RAM<double>(true);  
     testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}

