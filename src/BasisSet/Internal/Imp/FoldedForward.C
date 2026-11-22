// File: BasisSet/Internal/Imp/FoldedForward.C  The k-fold's multi-block forward: grouping by fold domain, and the folded call.
module;
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <array>     // foldDomain equality
module qchem.BasisSet.Internal.FoldedForward;
import qchem.BasisSet.PairProjection;   // PairProjection<dcmplx>: foldDomain, applyFolded, applyRawFolded

namespace qchem::BasisSet::Internal_Fold
{

namespace {

using Terms=std::vector<std::pair<const chmat_t*, const PairProjection<dcmplx>*>>;

bool CanFold(const PairProjection<dcmplx>& g, Route r)
{
    if (g.foldDomain[0]==0) return false;
    return r==Route::Ball ? bool(g.applyFolded) : bool(g.applyRawFolded);
}

//! A real TRIM block's D on the complex face the fold takes (upper triangle: hmat_t stores one).
chmat_t Widen(const hmat_t<double>& d)
{
    size_t n=d.rows();
    chmat_t c(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) c(i,j)=dcmplx(d(i,j),0.0);
    return c;
}

class FoldedForward : public MultiBlockForward
{
public:
    FoldedForward(Route r, size_t nBlocks, std::vector<const PairProjection<dcmplx>*> g, std::vector<std::vector<size_t>> groups,
                  std::vector<size_t> alone, std::unique_ptr<MultiBlockForward> aloneForward)
        : itsRoute(r), itsN(nBlocks), itsG(std::move(g)), itsGroups(std::move(groups)), itsAlone(std::move(alone))
        , itsAloneForward(std::move(aloneForward)) {}

    Projections Forward(const std::vector<FoldBlockD>& Ds) const override
    {
        if (Ds.size()!=itsN)
            throw std::logic_error("FoldedForward: "+std::to_string(Ds.size())+" density matrices for "+std::to_string(itsN)+" blocks");
        Projections ro;
        for (const auto& grp : itsGroups)
        {
            std::vector<chmat_t> widened; widened.reserve(grp.size());   // reserved: terms hold pointers into it
            Terms terms;
            for (size_t b : grp)
            {
                const FoldBlockD& d=Ds[b];
                if      (d.bloch) terms.push_back({d.bloch, itsG[b]});
                else if (d.real ) { widened.push_back(Widen(*d.real)); terms.push_back({&widened.back(), itsG[b]}); }
                else throw std::logic_error("FoldedForward: block "+std::to_string(b)+" has no density matrix");
            }
            const PairProjection<dcmplx>& lead=*itsG[grp.front()];
            if (itsRoute==Route::Ball) ro+=Projections(lead.applyFolded(terms));
            else                       ro+=Projections(lead.applyRawFolded(terms));
        }
        if (!itsAlone.empty())
        {
            std::vector<FoldBlockD> aloneD; aloneD.reserve(itsAlone.size());
            for (size_t b : itsAlone) aloneD.push_back(Ds[b]);
            ro+=itsAloneForward->Forward(aloneD);
        }
        return ro;
    }
private:
    Route                                      itsRoute;
    size_t                                     itsN;
    std::vector<const PairProjection<dcmplx>*> itsG;        //!< each block's overlap tensor (borrowed: the basis's integral cache)
    std::vector<std::vector<size_t>>           itsGroups;   //!< blocks of one fold domain, two or more
    std::vector<size_t>                        itsAlone;    //!< the rest, in order...
    std::unique_ptr<MultiBlockForward>         itsAloneForward;   //!< ...through the projector's per-block default
};

} //anon

std::unique_ptr<MultiBlockForward> MakeFoldedForward(const PairProjector<dcmplx>& p, const cFIT_SF_ABS& sf, Route route,
                                                     std::vector<FoldBlockRef> blocks)
{
    std::vector<const PairProjection<dcmplx>*> g(blocks.size(), nullptr);
    for (size_t b=0;b<blocks.size();b++)
    {
        if      (blocks[b].real ) g[b]=&blocks[b].real ->Overlap3C(sf);
        else if (blocks[b].bloch) g[b]=&blocks[b].bloch->Overlap3C(sf);
        else throw std::logic_error("MakeFoldedForward: block "+std::to_string(b)+" names no orbital block");
    }
    // Group by fold domain, in the blocks' order (a linear scan: a k mesh has few domains, and no map keyed on an identity).
    std::vector<std::vector<size_t>> groups;
    std::vector<size_t> alone;
    std::vector<bool> done(blocks.size(), false);
    for (size_t b=0;b<blocks.size();b++)
    {
        if (done[b]) continue;
        done[b]=true;
        if (!CanFold(*g[b], route)) { alone.push_back(b); continue; }
        std::vector<size_t> grp{b};
        for (size_t c=b+1;c<blocks.size();c++)
            if (!done[c] && g[c]->foldDomain==g[b]->foldDomain && CanFold(*g[c], route)) { grp.push_back(c); done[c]=true; }
        if (grp.size()>1) groups.push_back(std::move(grp)); else alone.push_back(b);
    }
    std::vector<FoldBlockRef> aloneBlocks; aloneBlocks.reserve(alone.size());
    for (size_t b : alone) aloneBlocks.push_back(blocks[b]);
    auto aloneForward=p.PairProjector<dcmplx>::ForwardBlocks(std::move(aloneBlocks));   // the base default, non-virtually
    return std::make_unique<FoldedForward>(route, blocks.size(), std::move(g), std::move(groups), std::move(alone), std::move(aloneForward));
}

} //namespace
