// File: Common/Imp/Parallel.C  The BLAS thread pin (rationale on the declaration).
module;
#include <cstdlib>   // setenv (the OpenMP wait policy -- see StopOmpThreadsBusyWaiting)
#include <cblas.h>   // openblas_set_num_threads -- an OpenBLAS extension, declared here by the
                     // openblas alternative of cblas.h (netlib's has no such call, which is the
                     // point: a BLAS swap breaks the build here rather than silently unpinning).
module qchem.Parallel;

namespace qchem {

void PinBlasToOneThread()
{
    openblas_set_num_threads(1);
}

void StopOmpThreadsBusyWaiting()
{
    // overwrite=0: an explicit setting in the environment WINS, so the A/B in the header stays runnable.
    setenv("KMP_BLOCKTIME",   "0",       0);   // libomp: spin 0 ms, then sleep (default is 200 ms)
    setenv("OMP_WAIT_POLICY", "PASSIVE", 0);   // the standard spelling, for any other runtime
}

} // namespace qchem
