// File: BasisSet/Imp/PairProjector.C  The periodic face's default multi-block forward (block by block, today's behaviour).
module;
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
module qchem.BasisSet.PairProjector;

namespace qchem::BasisSet
{

namespace {

//! The default \c ForwardBlocks: each block through its own \c Forward, summed (the additive identity adopts their realization).
class PerBlockForward : public MultiBlockForward
{
public:
    PerBlockForward(const PairProjector<dcmplx>& p, std::vector<FoldBlockRef> blocks) : itsP(p), itsBlocks(std::move(blocks)) {}
    Projections Forward(const std::vector<FoldBlockD>& Ds) const override
    {
        if (Ds.size()!=itsBlocks.size())
            throw std::logic_error("MultiBlockForward: "+std::to_string(Ds.size())+" density matrices for "+std::to_string(itsBlocks.size())+" blocks");
        Projections ro;
        for (size_t b=0;b<Ds.size();b++)
        {
            const FoldBlockRef& o=itsBlocks[b]; const FoldBlockD& d=Ds[b];
            if      (o.real  && d.real ) ro+=itsP.Forward(*o.real ).Forward(*d.real );
            else if (o.bloch && d.bloch) ro+=itsP.Forward(*o.bloch).Forward(*d.bloch);
            else throw std::logic_error("MultiBlockForward: block "+std::to_string(b)+"'s density matrix is not on the block's own face");
        }
        return ro;
    }
private:
    const PairProjector<dcmplx>& itsP;
    std::vector<FoldBlockRef>    itsBlocks;
};

} //anon

std::unique_ptr<MultiBlockForward> PairProjector<dcmplx>::ForwardBlocks(std::vector<FoldBlockRef> blocks) const
{
    return std::make_unique<PerBlockForward>(*this, std::move(blocks));
}

} //namespace
