// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"

// There is no process-wide integrals cache (doc/IntegralCacheOwnershipPlan.md S6): every cache has an owner, so nothing to set up here.

int main(int argc, char **argv)
{
     testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}


