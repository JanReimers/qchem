// File: BasisSet/Internal/FoldedForward.C  The k-FOLD's multi-block forward, shared by the projectors whose blocks can fold.
//
// One density's blocks on one lattice basis and one ladder (equal, non-zero PairProjection::foldDomain) collocate their SUMMED
// density in ONE box walk, through the first block's folded closure (user ruling 2026-11-22 "I agree A"; lane A's fold-k design
// note).  The G-space ρ̃ projector folds through \c applyFolded, the raw-raster projector through \c applyRawFolded; the grouping,
// the widening of a real TRIM block's D onto the complex face, and the per-block route for a block that cannot fold are the same,
// so they live here once.
module;
#include <memory>
#include <vector>
export module qchem.BasisSet.Internal.FoldedForward;
export import qchem.BasisSet.PairProjector;   // PairProjector<dcmplx>, FoldBlockRef, FoldBlockD, MultiBlockForward

export namespace qchem::BasisSet::Internal_Fold
{

//! Which folded closure of \c PairProjection<dcmplx> the forward runs: the G-ball (\c applyFolded, a ΔG_Map) or the raw raster
//! (\c applyRawFolded, an rvec_t).  The two projectors that fold differ only in this.
enum class Route {Ball, Raw};

//! \brief The folded forward for \a blocks of one density, as projector \a p (on fit basis \a sf) vends it: blocks of one fold
//! domain go through the \a route closure together; every other block goes through \a p's own per-block default
//! (\c PairProjector<dcmplx>::ForwardBlocks, called non-virtually), so the per-block loop exists once.
std::unique_ptr<MultiBlockForward> MakeFoldedForward(const PairProjector<dcmplx>& p, const cFIT_SF_ABS& sf, Route route,
                                                     std::vector<FoldBlockRef> blocks);

} //namespace
