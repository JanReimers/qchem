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
#include <string>    // ThreadSummary
export module qchem.Parallel;

export namespace qchem {

//! The pure parse behind every thread-count knob: a null or empty \a s gives \a dflt, anything else is atoi'd and
//! clamped to >=1 (so "0", "-3" and "abc" mean 1 -- serial -- rather than "all" or an error; the ruling on what 0
//! should mean is the open part of doc/CleanCode.md D-THREADS).  A function of its argument so it is unit-testable;
//! the env reads below happen once per process.
inline int ParseThreadCount(const char* s, int dflt)
{
    if (!s || !*s) return dflt;
    const int v=std::atoi(s);
    return v<1 ? 1 : v;
}

//! Worker threads for a parallel region: \c GPW_OMP_THREADS (read ONCE per process), clamped to >=1.
//! 1 (the default) means run serially -- callers keep a plain serial branch for it.
inline int WorkerThreads()
{
    static const int n=ParseThreadCount(std::getenv("GPW_OMP_THREADS"), 1);
    return n;
}

//! Thread CAP for the Becke-mesh build (\c UnitCell): \c GPW_OMP_THREADS when set, else 0 = "all the cores".
//! ⚠ This is the ONE region whose DEFAULT differs from \c WorkerThreads() (serial): its per-point partitions are
//! independent and slot-indexed, so the threaded build is bit-identical at any count and parallel-by-default costs
//! no anchor.  It is also why a run with \c GPW_OMP_THREADS unset can show ~500% CPU while \c WorkerThreads()==1
//! (found 2026-09-30) -- \c ThreadSummary() now states both on the run banner.
inline int MeshBuildThreads()
{
    // PRESERVES the historical reading: a value <1 ("0") here means ALL cores, whereas WorkerThreads() clamps it to 1.
    // That inconsistency is deliberate-until-ruled (doc/CleanCode.md D-THREADS: what should 0 mean?), not a design.
    static const int n=[]{ const char* s=std::getenv("GPW_OMP_THREADS"); const int v=s?std::atoi(s):0; return v<1 ? 0 : v; }();
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

//! The EFFECTIVE thread state of a run, one line, for the run banner: the pair/XC-loop workers, the Becke-mesh build,
//! and the BLAS count -- each with the knob that set it.  (A banner that printed only \c GPW_OMP_THREADS and a
//! hard-coded "BLAS pinned to 1" misdescribed a run whose mesh build used every core.)
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
