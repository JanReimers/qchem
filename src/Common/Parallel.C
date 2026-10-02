// File: Common/Parallel.C  The ONE opt-in worker-thread count, shared by every parallel region.
//
// The GPW pair loops (PG_Cart_MnD::NR_Evaluator::PairThreads) established the project's threading
// policy: SERIAL BY DEFAULT, opted into per run with QCHEM_OPENMP_THREADS (renamed 2026-10-02 from
// GPW_OMP_THREADS, which had nothing to do with basis sets; the old name is a deprecated alias).
// Serial-by-default is not
// timidity -- a threaded reduction sums in a load-dependent order, so the bit-anchors (and the
// OpenBLAS pin) only mean what they say on the serial path, and the suite runs many test binaries
// at once (ctest -j8) where a 16-thread fan-out per process would just thrash.
//
// This module is that knob, lifted out of the one evaluator that owned it, so the OTHER hot sites --
// the XC-mesh basis tables, the mesh-quadrature GEMMs, ... (doc/GPWPlan1.md item 1: "OMP coverage
// beyond the pair loops") -- read the SAME number instead of each growing its own getenv.  The name
// is QCHEM_OPENMP_THREADS.
//
// A caller that partitions by OUTPUT ELEMENT (each element still accumulated in one thread, in the
// serial order) is bit-identical at any thread count; one that partitions a REDUCTION is not, and
// must say so where it does it.
module;
#include <cstdlib>   // std::getenv/std::atoi
#include <string>    // ThreadSummary
export module qchem.Parallel;

export namespace qchem {

//! The number of PHYSICAL cores (not hardware threads): the distinct (package, core) pairs in sysfs
//! `thread_siblings_list`, else \c hardware_concurrency().  8 on an i7-10700 (16 hardware threads).  What
//! \c QCHEM_OPENMP_THREADS=0 ("auto") resolves to: SMT adds little for these floating-point loops, so half the
//! hardware threads is the right number only when SMT is on -- this reads the truth.
int PhysicalCores();

//! The pure parse behind the thread-count knob: null/empty \a s gives \a dflt; "0" means AUTO = \a autoCount;
//! a positive integer is itself; anything else (negative, garbage) is 1.  A function of its arguments so it is
//! unit-testable; the env read below happens once per process.
inline int ParseThreadCount(const char* s, int dflt, int autoCount)
{
    if (!s || !*s) return dflt;
    char* end=nullptr;
    const long v=std::strtol(s, &end, 10);
    if (end==s || *end!='\0' || v<0) return 1;
    return v==0 ? autoCount : int(v);
}

//! \brief THE thread count for every OpenMP region we own (pair loops, XC-mesh tables and quadrature, the
//! Becke-mesh build, the seed/FT sampling ...): \c QCHEM_OPENMP_THREADS, read ONCE per process.
//!   - unset  => 1: SERIAL.  Serial-by-default is not timidity: a threaded reduction sums in a load-dependent
//!     order, so the bit-anchors only mean what they say on the serial path, and the suite runs many test
//!     binaries at once (ctest -j8).
//!   - 0      => AUTO = \c PhysicalCores().  Opt-in only.
//!   - N>=1   => N.
//! ONE rule for every region (2026-10-02, user): the Becke-mesh build used to default to ALL cores while the others
//! defaulted to 1, which made a "serial" run show 500% CPU.  The old name \c GPW_OMP_THREADS (it had nothing to do
//! with basis sets) is still read as a DEPRECATED ALIAS when \c QCHEM_OPENMP_THREADS is absent.
inline int WorkerThreads()
{
    static const int n=[]
    {
        const char* s=std::getenv("QCHEM_OPENMP_THREADS");
        if (!s) s=std::getenv("GPW_OMP_THREADS");            // deprecated alias
        return ParseThreadCount(s, 1, PhysicalCores());
    }();
    return n;
}

//! \brief BLAS worker threads: \c QCHEM_BLAS_THREADS (read ONCE per process), clamped to >=1.
//! **Default 1** -- i.e. the historical pin, so every banked number is unchanged unless a run asks.
//!
//! ⚠ **DETERMINISM COMES FROM THE COUNT BEING FIXED, NOT FROM ITS BEING 1** (2026-09-06).  The
//! original pin was a correctness knob because OpenBLAS **auto-sizes its pool from machine load** when
//! left alone: the reduction order inside a GEMM then varies between runs of the same binary, and that
//! was measured moving an SCF total by >2e-5 with the machine-eps anchors flapping.  A FIXED count of
//! N is as deterministic as a fixed count of 1 -- it just sums in a different order, so it moves the
//! last ULP ONCE, as a re-bank, rather than run to run.
int BlasThreads();

//! The EFFECTIVE thread state of a run, one line, for the run banner: the OpenMP worker count and the BLAS count,
//! each with the knob that set it.
std::string ThreadSummary();

//! \brief Fix the BLAS to exactly \c BlasThreads() threads.  Call ONCE at the top of \c main() --
//! every test main and every CLI driver does.
//!
//! ⚠ **WHEN >1 IS SAFE, AND WHY IT IS NOT THE NESTING HAZARD IT LOOKS LIKE** (measured 2026-09-06).
//! The standing rule was "one level of parallelism, ours": flat OpenMP regions above (pair loops,
//! XC-mesh tables), BLAS pinned underneath.  But at the site that motivated this -- the XC-mesh
//! quadrature \f$H_{xc}=\Phi^\dagger\mathrm{diag}(wv)\Phi\f$ -- **there is no level above.**
//! `CompositeWF`'s Fock assembly walks the irrep blocks in a PLAIN SERIAL `for`, and there is no
//! `#pragma omp` anywhere in `WaveFunction` or `Hamiltonian`.  So the dispatched `zgemm` runs on one
//! core with eleven idle, and letting BLAS have them oversubscribes nothing: our OpenMP regions and
//! this GEMM never overlap in TIME.  (doc/Benchmark.md §7c read the same 1.21× as "the level above is
//! only 2 wide"; the level above is not 2 wide, it is not parallel at all.)
//! ⇒ Raise this only where that holds.  A future concurrent outer level over k-blocks/spins would make
//! `outer_width x BlasThreads() ~ cores` the rule instead, which is why the count is a NUMBER here and
//! not a bool.
//!
//! Deliberately a hard call into OpenBLAS rather than the \c OPENBLAS_NUM_THREADS env var (visible in
//! the source, not in someone's shell) and rather than a weak symbol (a BLAS swap must fail LOUDLY at
//! link time, not silently unpin the run).
void FixBlasThreads();

//! \brief Stop the OpenMP threads BUSY-WAITING between parallel regions.  Call ONCE at the top of
//! \c main(), beside \c FixBlasThreads -- same shape, same reason: a process-wide runtime setting
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
