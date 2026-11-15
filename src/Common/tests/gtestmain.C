// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"
import qchem.SCFCap;     // an SCF that ends at its iteration cap FAILS the running test (user ruling 2026-10-29)
#include <string>

// There is no process-wide integrals cache (doc/IntegralCacheOwnershipPlan.md S6): nothing to set up here.

int main(int argc, char **argv)
{
     testing::InitGoogleTest(&argc, argv);
    // EVERY SCF TEST FAILS IF AN SCF IN IT HITS NMaxIter UNCONVERGED (user 2026-10-29).  The SCF runs on the test's own
    // thread, so ADD_FAILURE lands in the running test.  A deliberately bounded loop says so with qchem::BoundedSCF.
    qchem::SetSCFCapHook([](const std::string& what){ ADD_FAILURE() << what; });
     return RUN_ALL_TESTS();
}


