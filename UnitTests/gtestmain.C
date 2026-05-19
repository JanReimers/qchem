// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"    
import qchem.BasisSet1.DB_Cache;

int main(int argc, char **argv)
{
      if (BasisSet1::theGlobalCache==0)
          BasisSet1::theGlobalCache=new BasisSet1::IntegralsCache_RAM<double>(true);  
     testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}

