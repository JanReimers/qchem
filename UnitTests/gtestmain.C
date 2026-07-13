// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"
#include <cblas.h>   // openblas_set_num_threads

// The integrals cache is now a construct-on-first-use singleton (BasisSet::theCache<T>()), so there is
// no longer any cache to new/delete here.

int main(int argc, char **argv)
{
    // REPRODUCIBILITY: pin OpenBLAS to one thread.  Left to itself OpenBLAS auto-sizes its internal thread pool
    // from machine load, and a multi-threaded BLAS/LAPACK reduction sums its terms in a load-dependent order --
    // so the last few ULPs of every eigenvalue/energy wobble from run to run (measured here: an SCF total energy
    // moving > 2e-5, and machine-eps regression anchors flapping between builds).  This is OpenBLAS's OWN internal
    // parallelism, not a thread-safety problem in our code: we call BLAS single-threaded, on private buffers.
    // Single-threading makes results bit-deterministic; the cost is tiny (our matrices are small, so ~4x of the
    // 5.2x GPW zgemm speed-up is SIMD, not threads).  Prefer this explicit call over an OPENBLAS_NUM_THREADS env
    // var so the choice is visible in the code, not hidden in a shell profile.
    openblas_set_num_threads(1);

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


