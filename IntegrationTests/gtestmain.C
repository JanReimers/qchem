// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"
import qchem.SCFCap;     // an SCF that ends at its iteration cap FAILS the running test (user ruling 2026-10-29)
#include <cstdlib>
#include <string>

// The integrals cache is now a construct-on-first-use singleton (BasisSet::theCache<T>()), so there is
// no longer any cache to new/delete here.

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    // EVERY SCF TEST FAILS IF AN SCF IN IT HITS NMaxIter UNCONVERGED (user 2026-10-29).  The SCF runs on the test's own
    // thread, so ADD_FAILURE lands in the running test.  A deliberately bounded loop says so with qchem::BoundedSCF.
    qchem::SetSCFCapHook([](const std::string& what){ ADD_FAILURE() << what; });
    // HERMETIC MATERIALS: a test sees the SHIPPED cells only, never the user's own library (~/.qchem/materials.json, OpenWork §2
    // 2026-11-09) -- a saved import must not move a name count or collide with a name a test states.  A test of that layer sets its own.
    setenv("QCHEM_USER_MATERIALS", "/nonexistent/qchem-test-user-materials.json", /*overwrite*/0);
    return RUN_ALL_TESTS();
}


