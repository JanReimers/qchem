// File: Imp/CompositeCD.C  Composite charge density: ONE set of irrep block densities over the full Irrep.
module;
#include <cassert>
#include <stdexcept>
#include <vector>
#include <memory>
#include <type_traits>
#include <variant>
#include <map>
#include <string>
module qchem.CompositeCD;
import qchem.ChargeDensity.Types;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO (each finite block's AO face)
import qchem.BasisSet.G_FieldEvaluator;  // G_RasterTransform -- the raster face that star-averages a grid field
import qchem.Blaze;
import qchem.Reporting;                 // report::Timed -- splitting the V_H field build (ParallelAndOraclePlan 1.3b)

namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> tComposite_CD<T>::tComposite_CD(std::vector<Symmetry::Lattice_3D::ReciprocalOp> pointOps)
    : itsPointOps(std::move(pointOps))
{};

// The VIEW ctor (V1.37): the parent's blocks, filtered by the caller, nothing owned.  Private -- only
// GetChannel builds one, and only over its own blocks, so a view can never outlive what it views.
template <class T> tComposite_CD<T>::tComposite_CD(std::vector<Block> blocks,
                                                   const std::vector<Symmetry::Lattice_3D::ReciprocalOp>& pointOps)
    : itsBlocks(std::move(blocks))
    , itsPointOps(pointOps)
{};

template <class T> tComposite_CD<T>::~tComposite_CD() {}

template <class T> void tComposite_CD<T>::Insert(std::unique_ptr<tDM_CD<double>> cd, const Irrep& qns)
{
    itsBlocks.push_back(Block{qns, cd_ref_t(cd.get())});
    itsOwned.push_back(cd_child_t(std::move(cd)));
    RebuildChannels();
}
template <class T> void tComposite_CD<T>::Insert(std::unique_ptr<tDM_CD<dcmplx>> cd, const Irrep& qns)
{
    itsBlocks.push_back(Block{qns, cd_ref_t(cd.get())});
    itsOwned.push_back(cd_child_t(std::move(cd)));
    RebuildChannels();
}

// The spin GROUPS of the walk: maximal runs of same-spin blocks, in walk order.  The wave function builds
// its blocks one spin irrep at a time (SpinIrreps(g) order: Up then Down, or just None), so a polarized
// composite has exactly two groups and a spin-agnostic one exactly one.  Only the HF sweep needs them
// (exchange is same-spin, and its Fock blocks are per SPATIAL irrep); everything else is a plain sum.
template <class T> std::vector<typename tComposite_CD<T>::Group> tComposite_CD<T>::Groups() const
{
    std::vector<Group> gs;
    for (size_t i=0;i<itsBlocks.size();++i)
        if (gs.empty() || gs.back().s!=itsBlocks[i].qns.ms) gs.push_back(Group{itsBlocks[i].qns.ms,i,i+1});
        else                                                gs.back().end=i+1;
    return gs;
}

// The channel views: one per polarized spin present, and only when there is something to filter OUT --
// a composite whose blocks all share one spin IS that channel and answers itself (see GetChannel), which
// is also what keeps a view from building views of its own.
template <class T> void tComposite_CD<T>::RebuildChannels()
{
    itsChannels.clear();
    const std::vector<Group> gs=Groups();
    if (gs.size()<2) return;
    for (const Group& g:gs)
    {
        if (!IsPolarized(g.s)) continue;
        std::vector<Block> mine;
        for (const Block& b:itsBlocks) if (b.qns.ms==g.s) mine.push_back(b);
        itsChannels[g.s].reset(new tComposite_CD(std::move(mine), itsPointOps));
    }
}

template <class T> const tChargeDensity<T>* tComposite_CD<T>::GetChannel(const Spin& s) const
{
    assert(IsPolarized(s) && "tComposite_CD::GetChannel: ask for Up or Down, not the total");
    const std::vector<Group> gs=Groups();
    if (gs.size()==1 && gs.front().s==s) return this;          // a single-spin set IS its channel
    auto i=itsChannels.find(s);
    return i==itsChannels.end() ? nullptr : i->second.get();   // null: this density does not resolve s
}

// Same-scalar view of a child slot, for the T-TYPED operations (contract clients, block maps, Phi tables):
// the argument's scalar is the composite FACE's T, so only a same-T child can consume it.  The cross arm is
// UNREACHABLE today (children are built by the same-T wave function); it becomes reachable -- and gets a
// genuine narrowing implementation instead of this throw -- when RealComplexPlan Step 3 gives blocks their
// own basis type.  A throw, not an assert: a wiring error here must be loud in Release too.
template <class T> static tDM_CD<T>& SameT(const cd_ref_t& c)
{
    if (auto* p=std::get_if<tDM_CD<T>*>(&c)) return **p;
    throw std::logic_error("tComposite_CD: a child block's scalar differs from the composite face -- "
                           "T-typed operations cannot cross scalars until RealComplexPlan Step 3 lands");
}

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//

// Whole-system Coulomb via ERI4 bra-ket symmetry (doc/ERI4Rework.md §4/§5.4).  Within ONE SPIN GROUP the
// composite holds one block per spatial irrep in the SAME order as Jall (both come from itsBS->Iterate, see
// CompositeWF::MakeIrrepWFs), so block k of the group <-> Jall[k].  ONE uniform loop over canonical pairs
// k<=l calls AccumulateDirectBoth, which handles the diagonal (k==l, self-pair -> single localized
// contraction) and the off-diagonal (ONE pass over the canonical J(k,l) scatters into BOTH Jall[k] and
// Jall[l], so J(l,k) is never fetched/built/cached) -- see IrrepCD::AccumulateDirectBoth.  Densities stay
// encapsulated in the IrrepCD leaves (the pair helper reaches its partner by a same-class cast, as
// MixIn/GetChangeFrom already do).  A POLARIZED composite runs the sweep once per spin group into the same
// Fock blocks -- Coulomb sees the TOTAL density; exchange is SAME-SPIN, so a pair never crosses a group.
// Cross-cast a composite leaf to the exact-exchange PAIR face (V1.6).  Every leaf on the real path is one;
// a failure means a density kind reached the HF sweep that is not a bra-ket pair partner, which used to be a
// silent no-op under -DNDEBUG (a zeroed J and a wrong Fock) and is now a loud error.
template <class T> static const tHF_Pair_CD<T>& PairOf(const tDM_CD<T>& cd)
{
    const tHF_Pair_CD<T>* p=dynamic_cast<const tHF_Pair_CD<T>*>(&cd);
    if (!p) throw std::runtime_error("Composite HF sweep: a density block is not a bra-ket pair partner "
                                     "(only an irrep LEAF is one) -- the whole-system J/K cannot be built.");
    return *p;
}

// V1.31: the WHOLE-SYSTEM route.  A basis with no per-irrep-pair ERI4 blocks (the SALC decorator -- R1.7)
// builds one whole-AO Fock and slices it, so driving it with the canonical-PAIR loop below made it rebuild
// the SAME matrix once per irrep.  J and K are LINEAR in D, so summing the blocks' AO densities FIRST and
// building ONCE is not an approximation -- it is the identity sum(F(D_C)) == F(sum(D_C)) -- and it turns N
// whole-AO builds into one.  Returns false when no leaf has the route (every ERI4 basis), leaving the pair
// loop untouched.  The probe is answered by the BASIS, once, at the top of the sweep: nothing per-pair, and
// nothing to memoize (this is what retired SymFockCache and its elementwise density compare).
// \a cds is ONE spin group's blocks (begin..end of the walk).
template <class T, class Blocks> static bool WholeSystemAll(const Blocks& cds, size_t begin, size_t end,
                                                            std::vector<hmat_t<T>>& Fall, bool exchange)
{
    if (begin==end) return false;
    const BasisSet::WholeSystemFock_IBS<T>* ws0=SameT<T>(cds[begin].cd).WholeSystemFock();
    if (!ws0) return false;
    hmat_t<T> Dao=blazem::zeroH<T>(ws0->AODimension());
    for (size_t k=begin;k<end;++k)
    {
        assert(SameT<T>(cds[k].cd).WholeSystemFock() && "mixed whole-system/pair bases in one composite density");
        SameT<T>(cds[k].cd).AddAODensity(Dao);
    }
    const hmat_t<T> Fao=ws0->MakeAOFock(Dao,exchange);      // the ONE build
    for (size_t k=begin;k<end;++k)
        SameT<T>(cds[k].cd).WholeSystemFock()->SliceAOFock(Fall[k-begin],Fao);
    return true;
}

// One spin group's sweep (Coulomb or exchange): the whole-system route if the basis has it, else the
// canonical-pair loop.  Fock block k of the group <-> Fall[k].
template <class Comp> void Composite_HFSystem<Comp>::SweepGroup(size_t begin, size_t end,
                                                                std::vector<hmat_t<double>>& Fall, bool exchange) const
{
    const auto& cds=self().itsBlocks;
    assert(Fall.size()==end-begin && "Fock blocks must be 1:1 with the spin group's irrep densities");
    if (WholeSystemAll<double>(cds,begin,end,Fall,exchange)) return;   // V1.31: one AO build + N slices, no pair loop
    for (size_t k=begin;k<end;++k)
        for (size_t l=k;l<end;++l)                                      // l>=k : diagonal + off-diagonal
        {
            const tHF_Pair_CD<double>& pk=PairOf(SameT<double>(cds[k].cd));
            const tHF_Pair_CD<double>& pl=PairOf(SameT<double>(cds[l].cd));
            if (exchange) pk.AccumulateExchangeBoth(Fall[k-begin],Fall[l-begin],pl);
            else          pk.AccumulateDirectBoth  (Fall[k-begin],Fall[l-begin],pl);
        }
}

template <class Comp> void Composite_HFSystem<Comp>::AccumulateDirectAll(std::vector<hmat_t<double>>& Jall) const
{
    for (const auto& g:self().Groups()) SweepGroup(g.begin, g.end, Jall, false);
}

// Exchange counterpart of AccumulateDirectAll (same canonical-pair structure; K(i,j)=K(j,i)^T).  Driven on
// ONE spin's channel by the polarized exchange term (exchange is same-spin), on the total by the RHF one --
// where the two groups of a polarized total sum into the same blocks (= K[D_total], then the term's -1/2).
template <class Comp> void Composite_HFSystem<Comp>::AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const
{
    for (const auto& g:self().Groups()) SweepGroup(g.begin, g.end, Kall, true);
}

template <class T> double tComposite_CD<T>::DM_ContractBlocks(const std::map<std::string,hmat_t<T>>& blocks) const
{
    double ret=0.0;
    for (const Block& b:itsBlocks) ret+=SameT<T>(b.cd).DM_ContractBlocks(blocks);
    return ret;
}

// rho on the quadrature = the sum of the blocks' contributions (each block's D already carries its BZ
// weight, so this is the k-average -- same convention as GetFourierDensity).  MIXED-RUN aware for free
// (3c-3): each child asks the quadrature with its OWN scalar and gets the matching table back, which is
// what replaced the two threaded table maps.  A cross-scalar (real) child inside a complex run is served
// by the same call it would make in a real run.
template <class T> rvec_t tComposite_CD<T>::ProjectOnto(const Fitting::ScalarProjector& p) const
{
    rvec_t ro(p.NumCoefficients(), 0.0);
    for (const Block& b:itsBlocks) std::visit([&](const auto* c){ ro+=c->ProjectOnto(p); }, b.cd);
    return ro;
}

// The energy contractions, MIXED-AWARE (Step 3c-2b): a same-scalar child contracts natively; a REAL
// child inside a complex run contracts through the term's real-block contract clients -- the static one
// IS a tStatic_CC<double> (same signature), the dynamic one is the run-typed Dynamic_CC_RealBlock the
// real periodic leaf consumes via its RealBlockEnergy_CD capability.
template <class T> double tComposite_CD<T>::DM_Contract(const tStatic_CC<T>* v) const
{
    double ret=0.0;
    for (const Block& blk:itsBlocks)
        ret+=std::visit([&](const auto* b)->double
        {
            using BT=std::remove_const_t<std::remove_pointer_t<std::decay_t<decltype(b)>>>;
            if constexpr (std::is_same_v<BT,tDM_CD<T>>) return b->DM_Contract(v);
            else if constexpr (std::is_same_v<T,dcmplx>)
            {
                auto* rc=dynamic_cast<const tStatic_CC<double>*>(v);
                if (!rc) throw std::logic_error("tComposite_CD: a real child needs the term's real-block "
                                                "contract client (Static_HT_RealBlock, RealComplexPlan 3c-1)");
                return b->DM_Contract(rc);
            }
            else throw std::logic_error("tComposite_CD: a complex child inside a real-faced run is impossible");
        }, blk.cd);
    return ret;
}

template <class T> double tComposite_CD<T>::DM_Contract(const tDynamic_CC<T>* v,const tDM_CD<T>* cd) const
{
    double ret=0.0;
    for (const Block& blk:itsBlocks)
        ret+=std::visit([&](const auto* b)->double
        {
            using BT=std::remove_const_t<std::remove_pointer_t<std::decay_t<decltype(b)>>>;
            if constexpr (std::is_same_v<BT,tDM_CD<T>>) return b->DM_Contract(v,cd);
            else if constexpr (std::is_same_v<T,dcmplx>)
            {
                auto* rcv=dynamic_cast<const Dynamic_CC_RealBlock*>(v);
                if (!rcv) throw std::logic_error("tComposite_CD: a real child needs the term's run-typed "
                                                 "energy client (Dynamic_CC_RealBlock, RealComplexPlan 3c-2b)");
                auto* rb=dynamic_cast<const RealBlockEnergy_CD*>(b);
                if (!rb) throw std::logic_error("tComposite_CD: a real child density must carry the "
                                                "RealBlockEnergy_CD capability (PeriodicIrrepCD<double>)");
                return rb->DM_ContractE(rcv, cd);
            }
            else throw std::logic_error("tComposite_CD: a complex child inside a real-faced run is impossible");
        }, blk.cd);
    return ret;
}

template <class T> double tComposite_CD<T>::GetTotalCharge() const
{
    double ret=0.0;
    for (const Block& b:itsBlocks) ret+=std::visit([](const auto* c){return c->GetTotalCharge();}, b.cd);
    return ret;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
// AO density-fit projection: sum the blocks' <rho|c>.  Each block is cross-cast to its AO face (finite
// path only -- a periodic composite is not a ProjectedDensity_AO, so the dcmplx body is inert), mirroring
// the FourierDensity cross-cast in GetFourierDensity below.
template <class T> rvec_t tComposite_CD<T>::GetRepulsion3C(const BasisSet::rFIT_CD_ABS* fbs) const
{
    if constexpr (std::is_same_v<T,double>)
    {
        rvec_t ret(fbs->GetNumFunctions(),0);
        for (const Block& b:itsBlocks)
            std::visit([&](const auto* c)
            {
                auto* ao=dynamic_cast<const Fitting::CoulombMetric_ProjectedDensity*>(c);
                assert(ao && "composite block has no Coulomb-metric projection face (finite path)");
                ret+=ao->GetRepulsion3C(fbs);
            }, b.cd);
        return ret;
    }
    else
        return rvec_t();   // inert: a periodic density carries no AO projection
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void tComposite_CD<T>::ReScale(double factor)
{
    // No UT coverage
    for (const Block& b:itsBlocks) std::visit([&](auto* c){c->ReScale(factor);}, b.cd);
    this->AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
}

// ★ WHY THIS THROWS (R2.5's remainder, closed 2026-09-09; generalized from the polarized container by
// V1.37).  `MixIn`/`GetChangeFrom` are BINARY operations on a hierarchy: `this` and `cd` must be the SAME
// representation -- a composite over the same blocks, in the same order -- before the algebra means
// anything, and no single-dispatch signature can express that.  (Double dispatch would move the same
// run-time check into a visitor and buy nothing.)  So the cast stays and it is the mixer's LINEAGE
// contract -- one mixer, one density family -- that the message names.  A throw rather than an assert
// because an assert is compiled out under NDEBUG, i.e. in every production run and benchmark -- exactly
// where a mixer fed the wrong lineage would actually happen.  Abstract->abstract, the intended idiom.
template <class T> static const tComposite_CD<T>& RequireCompositePartner(const tMixableDensity<T>& cd, const char* who)
{
    const tComposite_CD<T>* ecd=dynamic_cast<const tComposite_CD<T>*>(&cd);
    if (!ecd)
        throw std::runtime_error(std::string(who)+": the partner density is not a COMPOSITE density.  Mixing "
            "and convergence are BINARY operations -- both operands must be the same representation before "
            "the algebra means anything -- so a mixer holding a composite (per-irrep block) density can only "
            "be handed another one.  Reaching here means a mixer was seeded from one density family and "
            "driven with another, which is a composition error, not a recoverable condition.");
    return *ecd;
}

template <class T> void tComposite_CD<T>::MixIn(const tMixableDensity<T>& cd,double f)
{
    const tComposite_CD& ecd=RequireCompositePartner<T>(cd, "tComposite_CD::MixIn");
    if (itsBlocks.size()!=ecd.itsBlocks.size())
        throw std::runtime_error("tComposite_CD::MixIn: the two composites hold different block counts -- "
                                 "not the same irrep set (or the same imposed spin subgroup)");
    for (size_t i=0;i<itsBlocks.size();++i)
        std::visit([&](auto* mine, auto* theirs)
        {
            if constexpr (std::is_same_v<std::decay_t<decltype(mine)>,std::decay_t<decltype(theirs)>>)
                mine->MixIn(*theirs,f);
            else
                throw std::logic_error("tComposite_CD::MixIn: the two composites' child scalars differ per block");
        }, itsBlocks[i].cd, ecd.itsBlocks[i].cd);
    this->AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
}

template <class T> double tComposite_CD<T>::GetChangeFrom(const tMixableDensity<T>& cd) const
{
    const tComposite_CD& ecd=RequireCompositePartner<T>(cd, "tComposite_CD::GetChangeFrom");
    if (itsBlocks.size()!=ecd.itsBlocks.size())
        throw std::runtime_error("tComposite_CD::GetChangeFrom: the two composites hold different block counts "
                                 "-- not the same irrep set (or the same imposed spin subgroup)");
    double ret=0;
    for (size_t i=0;i<itsBlocks.size();++i)
        ret += std::visit([&](const auto* mine, const auto* theirs) -> double
        {
            if constexpr (std::is_same_v<std::decay_t<decltype(mine)>,std::decay_t<decltype(theirs)>>)
                return mine->GetChangeFrom(*theirs);
            else
                throw std::logic_error("tComposite_CD::GetChangeFrom: the two composites' child scalars differ per block");
        }, itsBlocks[i].cd, ecd.itsBlocks[i].cd);
    return ret;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double tComposite_CD<T>::operator()(const rvec3_t& r) const
{
    double ret=0.0;
    for (const Block& b:itsBlocks) ret+=std::visit([&](const auto* c){return c->operator()(r);}, b.cd);
    return ret;
}

template <class T> rvec3_t tComposite_CD<T>::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    rvec3_t ret(0,0,0);
    for (const Block& b:itsBlocks) ret+=std::visit([&](const auto* c){return c->Gradient(r);}, b.cd);
    return ret;
}

// The periodic face.  Each block's Fourier face, or a loud failure: every block on the plane-wave path is one.
static const FourierDensity& FourierOf(const cd_ref_t& c)
{
    const FourierDensity* fc=std::visit([](const auto* b){return dynamic_cast<const FourierDensity*>(b);}, c);
    assert(fc && "composite block is not a FourierDensity (plane-wave path)");
    return *fc;
}

// rho-tilde(Delta-m) = Sum_blocks rho-tilde_k (each block already BZ-weighted) = the BZ average Sum_k w_k rho_k.
template <class Comp> ΔG_Map Composite_Fourier<Comp>::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    ΔG_Map rg;
    for (const auto& b : self().itsBlocks)
        for (const auto& kv : FourierOf(b.cd).GetFourierDensity(c)) rg[kv.first]+=kv.second;
    return SymmetrizeGMap(rg, self().itsPointOps);   // IBZ star-average (no-op when {E}) -- doc/GPWPlan1.md item 3
}

// Raw rho_DM = Sum_k w_k rho_k(r) on c's ONE raster (the weights ride in each block's D).  ALL-OR-NOTHING:
// if any block lacks the raw path (returns empty) the whole composite answers empty, so the caller's E/H
// pair never mixes raw and ball blocks (doc/GPWPlan 0.5(f2)).
template <class Comp> rvec_t Composite_Fourier<Comp>::GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
{
    rvec_t sum;
    for (const auto& b : self().itsBlocks)
    {
        rvec_t r=FourierOf(b.cd).GetRhoOnGrid(c);
        if (r.size()==0) return rvec_t{};
        if (sum.size()==0) sum=std::move(r);
        else               sum+=r;
    }
    // IBZ: star-average the summed raster IN REAL SPACE (voxel permutation) so XC sees ρ_sym while staying
    // on the non-negative ρ_DM grid.  ASK THE RASTER (2026-08-24): the operation is a voxel permutation plus
    // an FFT glide, so it lives on G_RasterTransform, not on the fit face -- and this array IS a raster
    // array, so the cross-cast is the same "I want more" ask GetRhoOnGrid's whole route is built on.  A fit
    // basis without a raster cannot answer, and cannot have produced this array either.
    dynamic_cast<const BasisSet::G_RasterTransform&>(c).Symmetrize(sum);
    return sum;
}

// V_H = Sum_blocks V_H_k (V_H is linear in rho-tilde, so summing the per-block Coulomb projections == the
// projection of the summed density).  Each block bakes the kernel via its own Repulsion3C(c).
template <class Comp> ΔG_Map Composite_Fourier<Comp>::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    // ★ SPLIT INTO TWO BUCKETS (1.3b, 2026-09-06).  This whole call is 8.35 s of a 67 s MnO run and it
    // threads at 1.00×; doc/Benchmark.md §7c attributes it to the star-average, which had never been
    // measured apart from the per-block merge above it.  Both are ΔG_Map (std::map) walks, so both are
    // plausible -- hence two buckets rather than one guess.
    // ★ MERGE ALL BLOCKS RAW, THEN STAR-AVERAGE ONCE (2026-09-07, doc/ParallelAndOraclePlan.md 1.3b): a
    // polarized run used to pay the IBZ star-average TWICE for one field (once per channel).  The average
    // is linear, so Sym(up)+Sym(dn) is Sym(up+dn) exactly, and V_H depends on the TOTAL density anyway.
    // Measured on MnO: 6.5 s of a 63 s run at 0.042 s a call, so the duplicate was ~3 s of the wall.
    ΔG_Map rg=GetRepulsion3C_Raw(c);
    // V_H is linear in ρ̃ and |UG|=|G|, so symmetrizing V_H(G) == V_H of the symmetrized density -- exact.
    StarAverage(rg);
    return rg;
}

// The RAW half of the pair (FourierDensity): the per-block merge, WITHOUT the star-average.
template <class Comp> ΔG_Map Composite_Fourier<Comp>::GetRepulsion3C_Raw(const BasisSet::cFIT_CD_ABS& c) const
{
    qchem::report::Timed timed("scf: V_H per-block ΔG_Map merge");
    ΔG_Map rg;
    for (const auto& b : self().itsBlocks)
        for (const auto& kv : FourierOf(b.cd).GetRepulsion3C_Raw(c)) rg[kv.first]+=kv.second;
    return rg;
}

// The averaging half: THIS composite owns the crystal point ops, so it owns the star average.
template <class Comp> void Composite_Fourier<Comp>::StarAverage(ΔG_Map& rg) const
{
    if (self().itsPointOps.empty()) return;                 // {E}: exact no-op, and no bucket entry either
    qchem::report::Timed timed("scf: V_H IBZ star-average (SymmetrizeGMap)");
    rg=SymmetrizeGMap(rg, self().itsPointOps);
}

template class Composite_HFSystem<tComposite_CD<double>>;
template class tComposite_CD<double>;
template class Composite_Fourier<tComposite_CD<dcmplx>>;
template class tComposite_CD<dcmplx>;

} //namespace
