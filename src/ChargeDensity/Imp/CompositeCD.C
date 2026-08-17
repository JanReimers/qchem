// File: Imp/CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <cassert>
#include <stdexcept>
#include <vector>
#include <memory>
#include <type_traits>
#include <map>
#include <string>
module qchem.CompositeCD;
import qchem.ChargeDensity.Types;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO (each finite block's AO face)
import qchem.Blaze;

namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> tComposite_CD<T>::tComposite_CD(std::vector<Symmetry::Lattice_3D::ReciprocalOp> pointOps)
    : itsPointOps(std::move(pointOps))
{};

template <class T> void tComposite_CD<T>::Insert(std::unique_ptr<tDM_CD<T>> cd)
{
    itsCDs.push_back(std::move(cd));
}

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//

// Whole-system Coulomb via ERI4 bra-ket symmetry (doc/ERI4Rework.md §4/§5.4).  The composite holds one
// block per irrep in the SAME order as Jall (both come from itsBS->Iterate, see CompositeWF::MakeIrrepWFs),
// so block k <-> Jall[k].  ONE uniform loop over canonical pairs k<=l calls AccumulateDirectBoth, which
// handles the diagonal (k==l, self-pair -> single localized contraction) and the off-diagonal (ONE pass over
// the canonical J(k,l) scatters into BOTH Jall[k] and Jall[l], so J(l,k) is never fetched/built/cached) --
// see IrrepCD::AccumulateDirectBoth.  Densities stay encapsulated in the IrrepCD leaves (the pair helper
// reaches its partner by a same-class cast, as MixIn/GetChangeFrom already do).
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
template <class T, class CDV> static bool WholeSystemAll(const CDV& cds,
                                                         std::vector<hmat_t<T>>& Fall, bool exchange)
{
    if (cds.empty()) return false;
    const BasisSet::WholeSystemFock_IBS<T>* ws0=cds.front()->WholeSystemFock();
    if (!ws0) return false;
    hmat_t<T> Dao=blazem::zeroH<T>(ws0->AODimension());
    for (auto& c:cds)
    {
        assert(c->WholeSystemFock() && "mixed whole-system/pair bases in one composite density");
        c->AddAODensity(Dao);
    }
    const hmat_t<T> Fao=ws0->MakeAOFock(Dao,exchange);      // the ONE build
    for (size_t k=0;k<cds.size();++k)
        cds[k]->WholeSystemFock()->SliceAOFock(Fall[k],Fao);
    return true;
}

template <class Comp> void Composite_HFSystem<Comp>::AccumulateDirectAll(std::vector<hmat_t<double>>& Jall) const
{
    const auto& itsCDs=self().itsCDs;
    assert(Jall.size()==itsCDs.size() && "Fock blocks must be 1:1 with the composite irrep densities");
    if (WholeSystemAll<double>(itsCDs,Jall,false)) return;   // V1.31: one AO build + N slices, no pair loop
    const size_t N=itsCDs.size();
    for (size_t k=0;k<N;++k)
        for (size_t l=k;l<N;++l)                                            // l>=k : diagonal + off-diagonal
            PairOf(*itsCDs[k]).AccumulateDirectBoth(Jall[k],Jall[l],PairOf(*itsCDs[l]));
}

// Exchange counterpart of AccumulateDirectAll (same canonical-pair structure; K(i,j)=K(j,i)^T).  Driven
// per single-spin composite (exchange is same-spin), so the blocks it fills are one spin channel's K.
template <class Comp> void Composite_HFSystem<Comp>::AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const
{
    const auto& itsCDs=self().itsCDs;
    assert(Kall.size()==itsCDs.size() && "Fock blocks must be 1:1 with the composite irrep densities");
    if (WholeSystemAll<double>(itsCDs,Kall,true)) return;    // V1.31: exchange is linear in D too
    const size_t N=itsCDs.size();
    for (size_t k=0;k<N;++k)
        for (size_t l=k;l<N;++l)                                            // l>=k : diagonal + off-diagonal
            PairOf(*itsCDs[k]).AccumulateExchangeBoth(Kall[k],Kall[l],PairOf(*itsCDs[l]));
}

template <class T> double tComposite_CD<T>::DM_ContractBlocks(const std::map<std::string,hmat_t<T>>& blocks) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_ContractBlocks(blocks);
    return ret;
}

// rho at the caller's points = the sum of the blocks' contributions (each block's D already carries its
// BZ weight, so this is the k-average -- same convention as GetFourierDensity).
template <class T> rvec_t tComposite_CD<T>::DM_RhoAtPoints(const rvec3vec_t& r, const std::map<Irrep,mat_t<T>>& Phi) const
{
    rvec_t ro(r.size(), 0.0);
    for (auto& c:itsCDs) ro+=c->DM_RhoAtPoints(r,Phi);
    return ro;
}

template <class T> double tComposite_CD<T>::DM_Contract(const tStatic_CC<T>* v) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v);
    return ret;
}

template <class T> double tComposite_CD<T>::DM_Contract(const tDynamic_CC<T>* v,const tDM_CD<T>* cd) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v,cd);
    return ret;
}

template <class T> double tComposite_CD<T>::GetTotalCharge() const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->GetTotalCharge();
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
        for (auto& c:itsCDs)
        {
            auto* ao=dynamic_cast<const Fitting::CoulombMetric_ProjectedDensity*>(c.get());
            assert(ao && "composite block has no Coulomb-metric projection face (finite path)");
            ret+=ao->GetRepulsion3C(fbs);
        }
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
    for (auto& c:itsCDs) c->ReScale(factor);
    this->AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
}

template <class T> void tComposite_CD<T>::MixIn(const tMixableDensity<T>& cd,double f)
{
    const tComposite_CD* ecd = dynamic_cast<const tComposite_CD*>(&cd);
    assert(ecd);
    auto  b(ecd->itsCDs.begin());
    for (auto& c:itsCDs)
    {
        c->MixIn(**b,f);
        b++;
    }
    this->AdvanceHead();   // mutated in place -> Version() moved; keep this density the lineage head
}

template <class T> double tComposite_CD<T>::GetChangeFrom(const tMixableDensity<T>& cd) const
{
    const tComposite_CD* ecd = dynamic_cast<const tComposite_CD*>(&cd);
    assert(ecd);
    assert(itsCDs.size()==ecd->itsCDs.size());
    auto  b(ecd->itsCDs.begin());
    double ret=0;
    for (auto& c:itsCDs)
    {
        ret += c->GetChangeFrom(**b);
        b++;
    }
    return ret;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double tComposite_CD<T>::operator()(const rvec3_t& r) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->operator()(r);
    return ret;
}

template <class T> rvec3_t tComposite_CD<T>::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    rvec3_t ret(0,0,0);
    for (auto& c:itsCDs) ret+=c->Gradient(r);
    return ret;
}

// rho-tilde(Delta-m) = Sum_blocks rho-tilde_k (each block already BZ-weighted) = the BZ average Sum_k w_k rho_k.
template <class Comp> ΔG_Map Composite_Fourier<Comp>::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    ΔG_Map rg;
    for (const auto& blk : self().itsCDs)
    {
        auto* fc=dynamic_cast<const FourierDensity*>(blk.get());
        assert(fc && "composite block is not a FourierDensity (plane-wave path)");
        for (const auto& kv : fc->GetFourierDensity(c)) rg[kv.first]+=kv.second;
    }
    return SymmetrizeGMap(rg, self().itsPointOps);   // IBZ star-average (no-op when {E}) -- doc/GPWPlan1.md item 3
}

// Raw rho_DM = Sum_k w_k rho_k(r) on c's ONE raster (the weights ride in each block's D).  ALL-OR-NOTHING:
// if any block lacks the raw path (returns empty) the whole composite answers empty, so the caller's E/H
// pair never mixes raw and ball blocks (doc/GPWPlan 0.5(f2)).
template <class Comp> rvec_t Composite_Fourier<Comp>::GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
{
    rvec_t sum;
    for (const auto& blk : self().itsCDs)
    {
        auto* fc=dynamic_cast<const FourierDensity*>(blk.get());
        assert(fc && "composite block is not a FourierDensity (plane-wave path)");
        rvec_t r=fc->GetRhoOnGrid(c);
        if (r.size()==0) return rvec_t{};
        if (sum.size()==0) sum=std::move(r);
        else               sum+=r;
    }
    // IBZ: star-average the summed raster IN REAL SPACE (voxel permutation) so XC sees ρ_sym while staying
    // on the non-negative ρ_DM grid.  The fit basis owns its grid + the τ=0 direct ops; no-op unless folded.
    c.SymmetrizeRaster(sum);
    return sum;
}

// V_H = Sum_blocks V_H_k (V_H is linear in rho-tilde, so summing the per-block Coulomb projections == the
// projection of the summed density).  Each block bakes the kernel via its own Repulsion3C(c).
template <class Comp> ΔG_Map Composite_Fourier<Comp>::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    ΔG_Map rg;
    for (const auto& blk : self().itsCDs)
    {
        auto* fc=dynamic_cast<const FourierDensity*>(blk.get());
        assert(fc && "composite block is not a FourierDensity (plane-wave path)");
        for (const auto& kv : fc->GetRepulsion3C(c)) rg[kv.first]+=kv.second;
    }
    // V_H is linear in ρ̃ and |UG|=|G|, so symmetrizing V_H(G) == V_H of the symmetrized density -- exact.
    return SymmetrizeGMap(rg, self().itsPointOps);   // IBZ star-average (no-op when {E})
}

template class Composite_HFSystem<tComposite_CD<double>>;
template class tComposite_CD<double>;
template class Composite_Fourier<tComposite_CD<dcmplx>>;
template class tComposite_CD<dcmplx>;

} //namespace