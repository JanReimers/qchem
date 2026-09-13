// File: BasisSet/IntegralsCacheReport.C  The one thing an orchestrator may ask of the process-wide integral cache.
//
// V1.20d: the facades (Calculation, AtomCalculation) used to reach theCache<double>().EmitReport() -- the cache
// MECHANISM, an .Internal. module -- through a public face's re-export, i.e. a cross-family Internal import
// that the label was hiding.  What they want is exactly one sentence, so that sentence is the public face;
// the cache itself stays Internal.
module;
export module qchem.BasisSet.IntegralsCacheReport;

export namespace qchem::BasisSet
{
//! Snapshot the process-wide integrals cache (RAM tiers, hit/reuse rows) into the open run report's `cache`
//! section.  The cache is a singleton NEVER cleared between runs, so this is cumulative-to-this-point -- exact
//! for a one-run process, the running total in a multi-run one.  A no-op when no run is open.
//! \note An orchestrator-side call by nature: the cache's activity is continuous, and "what did the run
//! accumulate" only exists at the run's end, which the cache does not observe.
void EmitIntegralsCacheReport();
}
