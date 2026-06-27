// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"

// The integrals cache is now a construct-on-first-use singleton (BasisSet::theCache<T>()), so there is
// no longer any cache to new/delete here.

int main(int argc, char **argv)
{

     testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}


