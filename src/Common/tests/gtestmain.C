// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"
import qchem.Parallel;   // FixBlasThreads (the full rationale lives on the declaration)
import qchem.SCFCap;     // an SCF that ends at its iteration cap FAILS the running test (user ruling 2026-10-29)
#include <string>

// The integrals cache is now a construct-on-first-use singleton (BasisSet::theCache<T>()); nothing to
// new/delete here (this was a latent null-deref before -- the cache global was never created).

int main(int argc, char **argv)
{
     qchem::FixBlasThreads();
     qchem::StopOmpThreadsBusyWaiting();   // ONE level of parallelism (ours) + deterministic BLAS reductions
     testing::InitGoogleTest(&argc, argv);
    // EVERY SCF TEST FAILS IF AN SCF IN IT HITS NMaxIter UNCONVERGED (user 2026-10-29).  The SCF runs on the test's own
    // thread, so ADD_FAILURE lands in the running test.  A deliberately bounded loop says so with qchem::BoundedSCF.
    qchem::SetSCFCapHook([](const std::string& what){ ADD_FAILURE() << what; });
     return RUN_ALL_TESTS();
}


