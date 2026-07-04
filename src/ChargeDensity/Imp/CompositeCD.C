// File: Imp/CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <cassert>
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
template <class T> tComposite_CD<T>::tComposite_CD()
{};

template <class T> void tComposite_CD<T>::Insert(tDM_CD<T>* cd)
{
    itsCDs.push_back(std::unique_ptr<tDM_CD<T>>(cd));
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
template <class T> void tComposite_CD<T>::AccumulateDirectAll(std::vector<hmat_t<T>>& Jall) const
{
    assert(Jall.size()==itsCDs.size() && "Fock blocks must be 1:1 with the composite irrep densities");
    const size_t N=itsCDs.size();
    for (size_t k=0;k<N;++k)
        for (size_t l=k;l<N;++l)                                            // l>=k : diagonal + off-diagonal
            itsCDs[k]->AccumulateDirectBoth(Jall[k],Jall[l],*itsCDs[l]);
}

// Exchange counterpart of AccumulateDirectAll (same canonical-pair structure; K(i,j)=K(j,i)^T).  Driven
// per single-spin composite (exchange is same-spin), so the blocks it fills are one spin channel's K.
template <class T> void tComposite_CD<T>::AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall) const
{
    assert(Kall.size()==itsCDs.size() && "Fock blocks must be 1:1 with the composite irrep densities");
    const size_t N=itsCDs.size();
    for (size_t k=0;k<N;++k)
        for (size_t l=k;l<N;++l)                                            // l>=k : diagonal + off-diagonal
            itsCDs[k]->AccumulateExchangeBoth(Kall[k],Kall[l],*itsCDs[l]);
}

template <class T> double tComposite_CD<T>::DM_ContractBlocks(const std::map<std::string,hmat_t<T>>& blocks) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_ContractBlocks(blocks);
    return ret;
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
template <class T> rvec_t tComposite_CD<T>::GetRepulsion3C(const BasisSet::FIT_CD_ABS* fbs) const
{
    if constexpr (std::is_same_v<T,double>)
    {
        rvec_t ret(fbs->GetNumFunctions(),0);
        for (auto& c:itsCDs)
        {
            auto* ao=dynamic_cast<const Fitting::ProjectedDensity_AO*>(c.get());
            assert(ao && "composite block is not a ProjectedDensity_AO (finite path)");
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

template <class T> void tComposite_CD<T>::MixIn(const tDM_CD<T>& cd,double f)
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

template <class T> double tComposite_CD<T>::GetChangeFrom(const tDM_CD<T>& cd) const
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
template <class T> FourierMap tComposite_CD<T>::GetFourierDensity() const
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        FourierMap rg;
        for (const auto& c : itsCDs)
        {
            auto* fc=dynamic_cast<const FourierDensity*>(c.get());
            assert(fc && "composite block is not a FourierDensity (plane-wave path)");
            for (const auto& kv : fc->GetFourierDensity()) rg[kv.first]+=kv.second;
        }
        return rg;
    }
    else
    {
        assert(false && "a finite (non-periodic) density has no reciprocal-lattice Fourier series");
        return FourierMap{};
    }
}

template class tComposite_CD<double>;
template class tComposite_CD<dcmplx>;

} //namespace