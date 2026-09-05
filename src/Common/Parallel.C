// File: Common/Parallel.C  The ONE opt-in worker-thread count, shared by every parallel region.
//
// The GPW pair loops (PG_Cart_MnD::NR_Evaluator::PairThreads) established the project's threading
// policy: SERIAL BY DEFAULT, opted into per run with GPW_OMP_THREADS>1.  Serial-by-default is not
// timidity -- a threaded reduction sums in a load-dependent order, so the bit-anchors (and the
// OpenBLAS pin) only mean what they say on the serial path, and the suite runs many test binaries
// at once (ctest -j8) where a 16-thread fan-out per process would just thrash.
//
// This module is that knob, lifted out of the one evaluator that owned it, so the OTHER hot sites --
// the XC-mesh basis tables, the mesh-quadrature GEMMs, ... (doc/GPWPlan1.md item 1: "OMP coverage
// beyond the pair loops") -- read the SAME number instead of each growing its own getenv.  The name
// stays GPW_OMP_THREADS: it is what every production run script and plan doc already sets.
//
// A caller that partitions by OUTPUT ELEMENT (each element still accumulated in one thread, in the
// serial order) is bit-identical at any thread count; one that partitions a REDUCTION is not, and
// must say so where it does it.
module;
#include <cstdlib>   // std::getenv/std::atoi
export module qchem.Parallel;

export namespace qchem {

//! Worker threads for a parallel region: \c GPW_OMP_THREADS (read ONCE per process), clamped to >=1.
//! 1 (the default) means run serially -- callers keep a plain serial branch for it.
inline int WorkerThreads()
{
    static const int n=[]{ const char* s=std::getenv("GPW_OMP_THREADS"); const int v=s?std::atoi(s):1;
                           return v<1 ? 1 : v; }();
    return n;
}

//! \brief Pin the BLAS to ONE thread.  Call ONCE at the top of \c main() -- every test main and every
//! CLI driver does.
//!
//! The parallelism strategy is ONE level, ours: the flat OpenMP regions above (pair loops, XC-mesh
//! tables and GEMMs).  A BLAS that also threads sits UNDERNEATH those and does two damaging things:
//! it oversubscribes the box (our N workers x its M each), and -- the reason this is a correctness
//! knob, not a tuning one -- OpenBLAS auto-sizes its pool FROM MACHINE LOAD, so the reduction order
//! inside a GEMM becomes load-dependent and results drift in the last ULP between runs of the same
//! binary.  That was measured: an SCF total energy moved by >2e-5 and the machine-eps anchors
//! flapped.  Pinned, BLAS is deterministic and the threading we DO want stays where we can reason
//! about it.
//!
//! Deliberately a hard call into OpenBLAS rather than the \c OPENBLAS_NUM_THREADS env var (visible in
//! the source, not in someone's shell) and rather than a weak symbol (a BLAS swap must fail LOUDLY at
//! link time, not silently unpin the run).
void PinBlasToOneThread();

//! \brief Stop the OpenMP threads BUSY-WAITING between parallel regions.  Call ONCE at the top of
//! \c main(), beside \c PinBlasToOneThread -- same shape, same reason: a process-wide runtime setting
//! belongs in the source where it can be read, not in someone's shell.
//!
//! ⛔ WHY (measured 2026-09-04).  LLVM's libomp spins for **200 ms** after every parallel region before
//! letting a thread sleep (`KMP_BLOCKTIME`, default 200).  Our regions are PER SHELL PAIR -- short and
//! very numerous -- so the threads spend most of their life spinning on a barrier, and every spun cycle
//! is billed as CPU.  On NaF SR2 Γ at 12 threads:
//!
//!     default            wall 10.99 s   CPU 94.6 s
//!     KMP_BLOCKTIME=0    wall 10.92 s   CPU 32.6 s      <- same wall, same Etot, 62 s of pure spin gone
//!
//! ⇒ **65% of the billed CPU was doing nothing**, which made every threaded row in doc/Benchmark.md
//! meaningless -- the protocol's standing warning ("a 294 s serial build billed ~590 s CPU at 16
//! threads") was this, and §7's threaded table could not be filled while it stood.
//!
//! Set with overwrite=0, so an explicit `KMP_BLOCKTIME` / `OMP_WAIT_POLICY` in the environment still
//! wins -- the A/B above has to stay runnable.  Env rather than \c kmp_set_blocktime() because that is a
//! libomp extension needing <omp.h>, and this file deliberately makes no \c omp_*() calls.
//! \note MUST run before the first parallel region: libomp reads these at its own init, which is the
//! first OMP call in the process.  The top of \c main() is the only place that is guaranteed.
void StopOmpThreadsBusyWaiting();

} // namespace qchem
