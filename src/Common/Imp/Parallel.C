// File: Common/Imp/Parallel.C  The BLAS thread pin (rationale on the declaration).
module;
#include <fstream>
#include <set>
#include <string>
#include <thread>
#include <cstdlib>   // setenv (the OpenMP wait policy -- see StopOmpThreadsBusyWaiting)
#include <cblas.h>   // openblas_set_num_threads -- an OpenBLAS extension, declared here by the
                     // openblas alternative of cblas.h (netlib's has no such call, which is the
                     // point: a BLAS swap breaks the build here rather than silently unpinning).
module qchem.Parallel;

namespace qchem {

int BlasThreads()
{
    if (detail::g_blasOverride>0) return detail::g_blasOverride;      // a run's scope (ScopedWorkerThreads)
    static const int n=[]{ const char* s=std::getenv("QCHEM_BLAS_THREADS"); const int v=s?std::atoi(s):1;
                           return v<1 ? 1 : v; }();
    return n;
}

ScopedWorkerThreads::ScopedWorkerThreads(int openmp, int blas)
    : itsOld(detail::g_workerOverride), itsOldBlasOverride(detail::g_blasOverride), itsOldBlasPool(openblas_get_num_threads())
{
    if (openmp>0) detail::g_workerOverride=openmp;
    if (blas>0) { detail::g_blasOverride=blas; openblas_set_num_threads(blas); }
}

ScopedWorkerThreads::~ScopedWorkerThreads()
{
    if (detail::g_blasOverride!=itsOldBlasOverride || openblas_get_num_threads()!=itsOldBlasPool) openblas_set_num_threads(itsOldBlasPool);
    detail::g_blasOverride=itsOldBlasOverride;
    detail::g_workerOverride=itsOld;
}

int PhysicalCores()
{
    static const int n=[]
    {
        std::set<std::string> cores;
        for (int c=0; c<4096; ++c)
        {
            std::ifstream f("/sys/devices/system/cpu/cpu"+std::to_string(c)+"/topology/thread_siblings_list");
            if (!f) break;
            std::string sib; f>>sib;
            cores.insert(sib);                              // the SMT siblings of one core list the same string
        }
        const int hw=int(std::thread::hardware_concurrency());
        return cores.empty() ? (hw>0 ? hw : 1) : int(cores.size());
    }();
    return n;
}

std::string ThreadSummary()
{
    const char* q=std::getenv("QCHEM_OPENMP_THREADS");
    const char* g=std::getenv("GPW_OMP_THREADS");
    return "OpenMP regions=" + std::to_string(WorkerThreads()) + " ("
         + (detail::g_workerOverride>0 ? std::string("the run's threads.openmp")
            : q ? "QCHEM_OPENMP_THREADS="+std::string(q)
              : g ? "deprecated GPW_OMP_THREADS="+std::string(g)+" -- use QCHEM_OPENMP_THREADS"
                  : std::string("QCHEM_OPENMP_THREADS unset => serial")) + "), BLAS="
         // The count OpenBLAS ACTUALLY runs, not the request: an unpinned pool printed "BLAS=1" here while 16 threads ran (2026-10-31).
         + std::to_string(openblas_get_num_threads())
         + (detail::g_blasOverride>0 ? std::string(" (the run's threads.blas)")
            : " (QCHEM_BLAS_THREADS" + (std::getenv("QCHEM_BLAS_THREADS") ? "="+std::string(std::getenv("QCHEM_BLAS_THREADS")) : std::string(" unset => 1")) + ")")
         + (openblas_get_num_threads()!=BlasThreads() ? "  [WARNING: BLAS is NOT pinned to the requested "+std::to_string(BlasThreads())+"]" : "")
         // ⚠ N BLAS threads INSIDE each of N OpenMP threads is the nesting that measured as a net loss (oversubscription;
         // the OpenBLAS workers spin against our regions).  Say so on the line a benchmark row copies.
         + ((WorkerThreads()>1 && BlasThreads()>1) ? "  [WARNING: threaded BLAS under threaded OpenMP regions -- oversubscribes]" : "");
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

namespace {
// ★ THE PIN RUNS WHEN qcCommon LOADS, so EVERY executable gets it -- rundeck, gpwprobe, ITMain, the UT mains and the pybind
// module -- not only the mains that remember to call it.  rundeck and gpwprobe never did: OpenBLAS ran its default pool
// (16 spinning threads, ~1300% CPU beside a 6-thread run) while the banner said BLAS=1 (2026-10-31 perf of a β-MnO₂ PW run).
// Load order is safe: libopenblas is a dependency of qcCommon, so its own constructor (the pool) has run, and the libraries
// that call BLAS depend on qcCommon, so they load after this.  This is the ONE home: no main() calls these.
[[maybe_unused]] const bool theThreadPolicyIsSet=[]{ FixBlasThreads(); StopOmpThreadsBusyWaiting(); return true; }();
}

} // namespace qchem
