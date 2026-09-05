// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"
import qchem.Parallel;   // PinBlasToOneThread (the full rationale lives on the declaration)

// The integrals cache is now a construct-on-first-use singleton (BasisSet::theCache<T>()); nothing to
// new/delete here (this was a latent null-deref before -- the cache global was never created).

int main(int argc, char **argv)
{
     qchem::PinBlasToOneThread();
     qchem::StopOmpThreadsBusyWaiting();   // ONE level of parallelism (ours) + deterministic BLAS reductions
     testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
}


