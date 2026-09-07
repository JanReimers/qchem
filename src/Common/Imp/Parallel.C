// File: Common/Imp/Parallel.C  The BLAS thread pin (rationale on the declaration).
module;
#include <cstdlib>   // setenv (the OpenMP wait policy -- see StopOmpThreadsBusyWaiting)
#include <cblas.h>   // openblas_set_num_threads -- an OpenBLAS extension, declared here by the
                     // openblas alternative of cblas.h (netlib's has no such call, which is the
                     // point: a BLAS swap breaks the build here rather than silently unpinning).
module qchem.Parallel;

namespace qchem {

int BlasThreads()
{
    static const int n=[]{ const char* s=std::getenv("QCHEM_BLAS_THREADS"); const int v=s?std::atoi(s):1;
                           return v<1 ? 1 : v; }();
    return n;
}

void FixBlasThreads()
{
    // ⛔ STOP THE OpenBLAS POOL BUSY-WAITING FIRST -- it is the SAME defect StopOmpThreadsBusyWaiting
    // fixes for libomp, one runtime down, and without it a threaded BLAS is a NET LOSS (measured
    // 2026-09-06, MnO 12-thread, doc/ParallelAndOraclePlan.md 1.3):
    //
    //     QCHEM_BLAS_THREADS=1                        wall 67.07 s  CPU 587%  H_xc 9.10 s
    //     QCHEM_BLAS_THREADS=6                        wall 69.68 s  CPU 773%  H_xc 5.73 s   <- bucket
    //     QCHEM_BLAS_THREADS=12                       wall 76.97 s  CPU 905%  H_xc 6.11 s      1.6x
    //     QCHEM_BLAS_THREADS=6 + THREAD_TIMEOUT=1     wall 62.01 s  CPU 624%  H_xc 5.75 s      FASTER
    //
    // i.e. threading the GEMM made its own bucket 1.6x faster and the WHOLE RUN SLOWER, because the
    // OpenBLAS workers spin between calls and our OpenMP regions are what they steal from.  With the
    // spin off the bucket win reaches the wall (and V_H drops 8.35 -> 7.65 s for the same reason).
    // The value is an exponent-ish spin count, not ms; 1 is the practical minimum.
    //
    // MUST precede openblas_set_num_threads: OpenBLAS reads this when its pool initialises, and that
    // call is typically the first thing to initialise it.  overwrite=0, so an explicit setting in the
    // environment still wins and the A/B above stays runnable.
    setenv("OPENBLAS_THREAD_TIMEOUT", "1", 0);
    openblas_set_num_threads(BlasThreads());
}

void StopOmpThreadsBusyWaiting()
{
    // overwrite=0: an explicit setting in the environment WINS, so the A/B in the header stays runnable.
    setenv("KMP_BLOCKTIME",   "0",       0);   // libomp: spin 0 ms, then sleep (default is 200 ms)
    setenv("OMP_WAIT_POLICY", "PASSIVE", 0);   // the standard spelling, for any other runtime
}

} // namespace qchem
