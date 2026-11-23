// File: ChargeDensity/Internal/FoldableBlock.C  A density block that can join a k-FOLD: its orbital block and its D, for one call.
module;
#include <vector>
export module qchem.ChargeDensity.Internal.FoldableBlock;
export import qchem.BasisSet.PairProjector;   // FoldBlockRef, FoldBlockD, PairProjector<dcmplx>
export import qchem.ChargeDensity;            // DensityBlocks<T> (the whole-density pullback / virial)

export namespace qchem
{

//! \brief A periodic density BLOCK that can join its k-composite's fold (lane A 2026-11-22; user ruling "I agree A": the whole density
//! collocates itself).  The composite asks each of its blocks, through a cross-cast between abstract faces inside this library (the
//! leaves and the composite are one tightly coupled family, user 2026-11-22), and gathers the ones that answer into ONE
//! \c PairProjector<dcmplx>::ForwardBlocks call; a block that declines contracts by itself, as before.  Data-free.
//!
//! WHY THE LEAF DECIDES, NOT THE PROJECTOR: a leaf may project differently from the plain \f$D\f$ contraction (\c FactoredRho takes a
//! thin factor where the projector's forward \c FactorPays), and only the leaf knows; folding it anyway would change its projection.
class FoldableBlock
{
public:
    virtual ~FoldableBlock() = default;
    //! Append my orbital block to \a blocks and a reference to my D to \a Ds (valid while I am unchanged), and answer TRUE -- or answer
    //! FALSE, appending nothing, when projecting through \a p must stay my own call.
    virtual bool JoinFold(const BasisSet::PairProjector<dcmplx>& p, std::vector<BasisSet::FoldBlockRef>& blocks,
                          std::vector<BasisSet::FoldBlockD>& Ds) const=0;
    //! \brief Append my (block, spin, D) to \a b, for ONE whole-density pullback or virial (\c tStatic_PC::PullbackBlocks, user ruling
    //! 2026-11-22 "yes"), and answer TRUE -- or FALSE, appending nothing, when my derivative must stay my own \c DM_Pullback call.  A real
    //! block in a complex run joins \a b's \c real list.
    virtual bool JoinPullback(ChargeDensity::DensityBlocks<double>& b) const=0;
    virtual bool JoinPullback(ChargeDensity::DensityBlocks<dcmplx>& b) const=0;
};

} // namespace qchem
