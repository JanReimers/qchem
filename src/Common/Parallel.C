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

} // namespace qchem
