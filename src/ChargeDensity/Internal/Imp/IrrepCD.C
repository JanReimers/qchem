// File: ExactIrrepCD.C  Exact implementation of the charged density.
module;
#include <algorithm>   // std::min (the threaded rho-sampling block partition)
#include <cassert>
#include <complex>
#include <exception>   // std::exception_ptr (throw containment across the threaded block loop)
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <type_traits>
#include <vector>
#include <atomic>
#include <map>
#include <string>

module qchem.ChargeDensity.Imp.IrrepCD;
import qchem.Symmetry;
import qchem.Blaze;
import qchem.Parallel;                  // WorkerThreads (GPW_OMP_THREADS -- the rho-sampling GEMM)
import qchem.BasisSet.Orbital_DFT_IBS;   // cast the basis UP to the G-space capability (dcmplx path)

namespace qchem::ChargeDensity
{

// The density-freshness serial source (Version()) now lives in qchem.ChargeDensity as a SINGLE shared
// counter (NextDensityVersion): every density kind must draw from the same global clock or their serials
// collide across kinds and a dynamic term reuses a stale cached matrix (see the module doc + the SAD seed
// types NumericCD/SeedCD, which also stamp from it).

typedef Vector3D<std::complex<double> > Vec3;

rvec3_t  GradientContraction(const vec_t<rvec3_t >&, const vec_t<double>&, const rsmat_t&);
// rvec3_t  GradientContraction(const Vector<Vec3 >&, const Vector<std::complex<double> >&, const smat_t<std::complex<double> >&);

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> IrrepCD<T>::IrrepCD()
    : itsVersion(NextDensityVersion())
{};

template <class T> IrrepCD<T>::IrrepCD(const DenSMat& D,const tobs_t<T>* theBasisSet,Irrep qns)
    : itsDensityMatrix(D)
    , itsBasisSet(theBasisSet)
    , itsSpin(qns.ms)
    , itsIrrep(qns)
    , itsVersion(NextDensityVersion())
{
    assert(itsBasisSet);
};

template <class T> bool IrrepCD<T>::IsZero() const
{
    // return max(abs(itsDensityMatrix))==0.0;
    return blazem::isZero(itsDensityMatrix);
}


//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
template <> void IrrepCD<double>::AccumulateDirect(rsmat_t& Jii) const
{
    const rohfbs_t* bs=dynamic_cast<const rohfbs_t*>(itsBasisSet);
    assert(bs);
    if (!IsZero()) bs->AccumulateDirect(Jii,itsDensityMatrix,bs);   // diagonal: this block's basis on both sides
}

template <> void IrrepCD<double>::AccumulateExchange(rsmat_t& Kii) const
{
    const rohfbs_t* bs=dynamic_cast<const rohfbs_t*>(itsBasisSet);
    assert(bs);
    if (!IsZero()) bs->AccumulateExchange(Kii,itsDensityMatrix,bs);
}

// One canonical irrep pair (this=i, other=j, i<=j) scattered into both irreps' Coulomb blocks.  This is the
// SOLE entry point the composite loops over (i<=j, so the DIAGONAL i==j lands here too): when other IS this
// (the self-pair) there is no bra-ket partner, so it collapses to a single localized contraction --
// ScatterBoth on the diagonal would add J.D + J^T.D = 2 J.D (the block is bra-ket symmetric).  Off-diagonal:
// the partner density is reached by a same-class cast (the IrrepCD<->IrrepCD idiom used by MixIn /
// GetChangeFrom); both empty -> nothing to build; the basis fetches ONLY the canonical J(i,j) block, so
// J(j,i) is never materialized.
template <class Leaf> void IrrepCD_HFPair<Leaf>::AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const tHF_Pair_CD<double>& other) const
{
    if (&other==this) { self().AccumulateDirect(Ji); return; }  // diagonal self-pair: single-block contraction
    // V1.8: no cast to the partner's concrete type.  We hand IT our block and let it finish with its own --
    // the contraction belongs to the BASIS, and a partner's only job is to present its matrix to that
    // routine.  Asking "give me your density matrix" would have reintroduced the getter IrrepCD is
    // deliberately without (CLAUDE.md cites this very class).
    const rohfbs_t* bs_i=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_i && "IrrepCD::AccumulateDirectBoth: this block's basis has no HF contraction face");
    other.CompleteDirectPair(Ji,Jj,self().itsDensityMatrix,bs_i);
}

// The partner's half: our own block completes the caller's pair.  Ji/Jj and Di/bs_i belong to the CALLER
// (irrep i); this object is irrep j.  The zero pre-screen lives here because only now are both blocks known.
template <class Leaf> void IrrepCD_HFPair<Leaf>::CompleteDirectPair(rsmat_t& Ji, rsmat_t& Jj, const rsmat_t& Di,
                                                    const BasisSet::Orbital_HF_IBS<double>* bs_i) const
{
    if (blazem::isZero(Di) && self().IsZero()) return;          // both empty: nothing to build
    const rohfbs_t* bs_j=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_j && "IrrepCD::CompleteDirectPair: this block's basis has no HF contraction face");
    bs_i->AccumulateDirectBoth(Ji,Jj,Di,self().itsDensityMatrix,bs_j);
}

// --- V1.31 whole-system route ---------------------------------------------------------------------------
// The basis answers whether it has one; this block only has to fold its own density up to the AO space.  The
// SLICE back down needs no density, so the composite drives that through the basis face directly.
template <> const BasisSet::WholeSystemFock_IBS<double>* IrrepCD<double>::WholeSystemFock() const
{
    return dynamic_cast<const BasisSet::WholeSystemFock_IBS<double>*>(itsBasisSet);
}

template <> void IrrepCD<double>::AddAODensity(rsmat_t& Dao) const
{
    const BasisSet::WholeSystemFock_IBS<double>* ws=WholeSystemFock();
    assert(ws && "IrrepCD::AddAODensity: this block's basis has no whole-system Fock route");
    if (!IsZero()) ws->AddAODensity(Dao,itsDensityMatrix);
}

// Exchange counterpart (see AccumulateDirectBoth): diagonal -> single AccumulateExchange, else the canonical
// Exchange block is fetched once.
template <class Leaf> void IrrepCD_HFPair<Leaf>::AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const tHF_Pair_CD<double>& other) const
{
    if (&other==this) { self().AccumulateExchange(Ki); return; }  // diagonal self-pair: single-block contraction
    const rohfbs_t* bs_i=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_i && "IrrepCD::AccumulateExchangeBoth: this block's basis has no HF contraction face");
    other.CompleteExchangePair(Ki,Kj,self().itsDensityMatrix,bs_i);
}

template <class Leaf> void IrrepCD_HFPair<Leaf>::CompleteExchangePair(rsmat_t& Ki, rsmat_t& Kj, const rsmat_t& Di,
                                                      const BasisSet::Orbital_HF_IBS<double>* bs_i) const
{
    if (blazem::isZero(Di) && self().IsZero()) return;
    const rohfbs_t* bs_j=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_j && "IrrepCD::CompleteExchangePair: this block's basis has no HF contraction face");
    bs_i->AccumulateExchangeBoth(Ki,Kj,Di,self().itsDensityMatrix,bs_j);
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
// AO density-fit projection <rho|c> = Sum_ab D_ab <ab|c>, the finite (double) path's ProjectedDensity_AO
// face.  The periodic (dcmplx) density is NOT a ProjectedDensity_AO (see ProjectedDensityBase), so this is
// never reached for dcmplx; the if-constexpr keeps the double-only 3-centre machinery out of that build.
template <class T> rvec_t IrrepCD<T>::GetRepulsion3C(const BasisSet::rFIT_CD_ABS* fbs) const
{
    if constexpr (std::is_same_v<T,double>)
    {
        if (IsZero()) return rvec_t(fbs->GetNumFunctions(),0.0);
        auto dftbs=dynamic_cast<const todftbs_t<T>*>(itsBasisSet);
        assert(dftbs);
        // Contract the density matrix against the basis's CACHED, D-free 3-centre projection tensor <ab|c>
        // HERE -- the DENSITY owns D, so the D-contraction is a density operation, not a basis one.  The
        // basis exposes only the tensor (Repulsion3C(c) -> Projector3, built once, keyed by BasisSetID);
        // D never crosses into qcBasisSet.  This is the real-space model for fixing MakeFourierDensity(D): the
        // {G} 3-centre integral is the delta <ij|Dm>, and rho-tilde = Sum_ij D_ij <ij|Dm> is the SAME contraction.
        const auto& R=dftbs->Repulsion3C(*fbs).dense;     // <ab|c> (dense realization: one smat per fit function c)
        rvec_t ret(fbs->GetNumFunctions());
        for (size_t i=0;i<R.size();++i)
            ret[i]=blazem::sum(itsDensityMatrix % R[i]);  // <rho|c_i> = Sum_ab D_ab <ab|c_i>
        return ret;
    }
    else
        return rvec_t();   // inert: a periodic density carries no AO projection
}


// Energy = trace(D V) = Sum_ij D_ij V_ji = Sum_ij D_ij trans(V)_ij = sum(D % trans(V)).  For real
// symmetric V this is identical to sum(D % V) (the old form); for complex HERMITIAN D,V the transpose
// matters -- sum(D % V) would silently give the wrong (still-real) value.
template <class T> double IrrepCD<T>::DM_Contract(const tStatic_CC<T>* v) const
{
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(v->GetMatrix(itsBasisSet,itsSpin)));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

template <class T> double IrrepCD<T>::DM_Contract(const tDynamic_CC<T>* v,const tDM_CD<T>* cd) const
{
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(v->GetEMatrix(itsBasisSet,itsSpin,cd)));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

// This irrep's contribution to a whole-system energy: sum(D % B_i^T) with B_i the block for this basis.
// Same contraction as DM_Contract above, but the Fock block comes from a caller-supplied map (the HF
// term's cached blocks) instead of a per-irrep GetMatrix call -- so no GetMatrix round-trip.
template <class T> double IrrepCD<T>::DM_ContractBlocks(const std::map<std::string,hmat_t<T>>& blocks) const
{
    if (IsZero()) return 0.0;
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(blocks.at(itsBasisSet->BasisSetID())));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

// This irrep's rho at the caller's points from its pre-evaluated basis table (the FIELD dual of
// DM_ContractBlocks above): rho_g = Re[Phi D Phi^dag]_gg = Re sum_j (Phi D)_gj conj(Phi_gj) -- one
// (npts x n)(n x n) GEMM + a cheap rowwise dot, matching operator()(r)'s trans(phi) D conj(phi) per row.
// No table for this basis -> self-evaluate pointwise (correct; the caller's first pass may not have
// built every block's table yet).
template <class T> rvec_t IrrepCD<T>::DM_RhoAtPoints(const rvec3vec_t& r, const std::map<Irrep,mat_t<T>>& Phi) const
{
    rvec_t ro(r.size(), 0.0);
    if (IsZero()) return ro;
    // Spatial key: Phi is a property of the basis block alone, so the two spin channels of one
    // spatial block share one table (Irrep interleaves spin into its ordering, hence Spin::None).
    auto it=Phi.find(Irrep(Spin::None,itsIrrep.sym));
    if (it==Phi.end())
    {
        for (size_t g=0; g<r.size(); g++) ro[g]=(*this)(r[g]);
        return ro;
    }
    const mat_t<T>& P=it->second;
    assert(P.rows()==r.size());
    assert(P.columns()==itsDensityMatrix.rows());
    // rho_g = Phi_g^dag D Phi_g, as (P D) row-dotted back into P: O(npts n^2), paid once per density
    // serial and THE largest per-iteration bucket of the atom-centred XC route (MnO Γ: 12 s/iteration
    // at 48k points x 122 functions).  Threaded by MESH-POINT BLOCK: rho_g depends on row g alone, so
    // each block is an independent GEMM + dot and the result is bit-identical at any thread count (a
    // partition of the OUTPUT, not of a reduction).  Opt-in via GPW_OMP_THREADS (qchem.Parallel).
    // D as a PLAIN matrix: an adaptor operand (HermitianMatrix) is not BLAS-dispatchable, so under
    // BLAZE_BLAS_MODE the product would silently fall back to blaze's own kernel -- 16x slower than
    // zgemm on this shape (measured).  One n x n copy per call is nothing beside the npts x n x n GEMM.
    const mat_t<T> D(itsDensityMatrix);
    auto block=[&](size_t g0, size_t g1)
    {
        auto Pb = blazem::submatrix(P, g0, size_t(0), g1-g0, P.columns());
        mat_t<T> PD = Pb*D;
        for (size_t g=g0; g<g1; g++)
        {
            T acc{};
            for (size_t j=0; j<PD.columns(); j++)
                if constexpr (std::is_same_v<T,dcmplx>) acc+=PD(g-g0,j)*std::conj(P(g,j));
                else                                    acc+=PD(g-g0,j)*P(g,j);
            ro[g]=std::real(acc);
        }
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1 && r.size()>size_t(nthreads))
    {
        const size_t blk=(r.size()+size_t(nthreads)-1)/size_t(nthreads);
        std::exception_ptr firstEx;                   // throw containment (an escape = std::terminate)
        #pragma omp parallel for schedule(static,1) num_threads(nthreads)
        for (size_t b=0; b<(r.size()+blk-1)/blk; b++)
        {
            try { block(b*blk, std::min(r.size(), (b+1)*blk)); }
            catch (...)
            {
                #pragma omp critical (cd_rho_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return ro;
    }
#endif
    // Serial: multiply the WHOLE table, not a full-height view -- a submatrix operand costs blaze its
    // aligned/padded kernel choice, and this is the default (unthreaded) path every suite run takes.
    mat_t<T> PD = P*D;
    for (size_t g=0; g<r.size(); g++)
    {
        T acc{};
        for (size_t j=0; j<PD.columns(); j++)
            if constexpr (std::is_same_v<T,dcmplx>) acc+=PD(g,j)*std::conj(P(g,j));
            else                                    acc+=PD(g,j)*P(g,j);
        ro[g]=std::real(acc);
    }
    return ro;
}

template <class T> double IrrepCD<T>::GetTotalCharge() const
{
    // N = integral rho = Sum_ij D_ij S_ji = Tr(D S) = sum(D % trans(S)) (% is the blaze Schur/direct product).
    // The trans matters for a genuinely COMPLEX Hermitian D,S (the complex-k Bloch case): sum(D % S) would be
    // Tr(D S^T) -- a different (still-real) value.  Same transpose fix as DM_Contract above; a no-op for real
    // symmetric S (trans(S)==S, byte-identical at Gamma / half-integer k).
    return std::real(blazem::sum(itsDensityMatrix%blazem::trans(itsBasisSet->Overlap())));
}


//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void IrrepCD<T>::ReScale(double factor)
{
    // No UT coverage
    itsDensityMatrix*=factor;
    itsVersion=NextDensityVersion();   // a mutation is logically a new density
}

template <class T> void IrrepCD<T>::MixIn(const tMixableDensity<T>& cd,double c)
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    itsDensityMatrix = itsDensityMatrix*(1-c) + eicd->itsDensityMatrix*c;
    itsVersion=NextDensityVersion();   // a mutation is logically a new density
}

template <class T> double IrrepCD<T>::GetChangeFrom(const tMixableDensity<T>& cd) const
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    return std::real(blazem::norm(itsDensityMatrix - eicd->itsDensityMatrix));
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double IrrepCD<T>::operator()(const rvec3_t& r) const
{
    vec_t<T> phir=(*itsBasisSet)(r);
    return std::real(blazem::trans(phir)*itsDensityMatrix*blazem::conj(phir));
}

template <class T> rvec3_t IrrepCD<T>::Gradient(const rvec3_t& r) const
{
    // No UT coverage
    vec_t<T> phir=(*itsBasisSet)(r);
    vec_t<rvec3_t > gphir=itsBasisSet->Gradient(r);
    return GradientContraction(gphir,phir,itsDensityMatrix);
}

// rho-tilde from the density matrix, via the basis's G-space capability (plane-wave / dcmplx only).
// itsDensityMatrix already carries the BZ weight w_k (TOrbitals::GetChargeDensity scales it), so the
// composite's sum over blocks is the BZ average Sum_k w_k rho_k.
template <class Leaf> ΔG_Map IrrepCD_Fourier<Leaf>::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    // Contract D against the basis's D-free OVERLAP tensor (empty kernel) HERE -- the DENSITY owns D, so
    // this is where rho-tilde(dm) = (1/Omega) Sum_{G_i-G_j=dm} D_ij is formed.  The overlap-metric sibling
    // of GetRepulsion3C above; D never crosses into the basis.
    auto* fb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx>*>(self().itsBasisSet);
    assert(fb && "GetFourierDensity requires a Orbital_DFT_IBS<dcmplx> (plane-wave) basis");
    return Contract(fb->Overlap3C(c), self().itsDensityMatrix);
}

// Raw rho_DM(r) on c's raster (doc/GPWPlan 0.5(f2)): the basis's collocation-native forward, when it has
// one (GPW's Overlap3C carries applyRaw; a plane-wave basis leaves it empty -> we answer empty and the
// caller falls back to the ball route).  itsDensityMatrix carries the BZ weight, exactly as above.
template <class Leaf> rvec_t IrrepCD_Fourier<Leaf>::GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
{
    auto* fb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx>*>(self().itsBasisSet);
    assert(fb && "GetRhoOnGrid requires a Orbital_DFT_IBS<dcmplx> (plane-wave) basis");
    const Projector3<dcmplx>& g=fb->Overlap3C(c);
    return g.applyRaw ? g.applyRaw(self().itsDensityMatrix) : rvec_t{};
}

// V_H(dm) = 4pi rho-tilde(dm)/|G|^2: contract D against the basis's D-free COULOMB tensor (Repulsion3C(c), the
// diagonal kernel baked in) -- the reciprocal mirror of the finite GetRepulsion3C(fbs) above.  D stays here.
template <class Leaf> ΔG_Map IrrepCD_Fourier<Leaf>::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    auto* fb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx>*>(self().itsBasisSet);
    assert(fb && "GetRepulsion3C requires a Orbital_DFT_IBS<dcmplx> (plane-wave) basis");
    return Contract(fb->Repulsion3C(c), self().itsDensityMatrix);
}


//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& IrrepCD<T>::Write(std::ostream& os) const
{
    return os << itsDensityMatrix;
}

template class IrrepCD_HFPair<IrrepCD<double>>;
template class IrrepCD<double>;

// --- Complex (plane-wave) density.  HF does NOT apply (the plane-wave path uses no 4-centre exchange), so
// those are NA; the gradient (GGA/plotting) is not yet wired for complex and is unused by the LDA SCF.  The
// AO density-fit projection (GetRepulsion3C) is no longer forced on this path -- IrrepCD<dcmplx> is not a
// ProjectedDensity_AO -- so its old NA-assert specialization is gone.  The remaining members are the generic
// templates above (DM_Contract, op()(r), GetTotalCharge, MixIn, GetChangeFrom, ...), complex-correct.
template <> void IrrepCD<dcmplx>::AccumulateDirect(hmat_t<dcmplx>&) const
{ assert(false && "AccumulateDirect: HF not applicable to a complex plane-wave density"); }
template <> void IrrepCD<dcmplx>::AccumulateExchange(hmat_t<dcmplx>&) const
{ assert(false && "AccumulateExchange: HF not applicable to a complex plane-wave density"); }
template <> const BasisSet::WholeSystemFock_IBS<dcmplx>* IrrepCD<dcmplx>::WholeSystemFock() const
{ return nullptr; }   // no SALC decorator on the periodic path: the pair route is the only one
template <> void IrrepCD<dcmplx>::AddAODensity(hmat_t<dcmplx>&) const
{ assert(false && "AddAODensity: HF not applicable to a complex plane-wave density"); }
// NOT IMPLEMENTED (not "zero"): the complex GradientContraction the real path uses has no dcmplx sibling,
// because nothing on the periodic path is a GGA or a gradient plotter.  Throwing keeps the FIRST such
// consumer honest -- the silent rvec3_t(0,0,0) this replaced handed it a plausible wrong field (R1.4).
template <> rvec3_t IrrepCD<dcmplx>::Gradient(const rvec3_t&) const
{
    throw std::logic_error("IrrepCD<dcmplx>::Gradient: grad(rho) is not wired for the complex (periodic) "
                           "density -- the periodic path is LDA-only.  Implement the complex gradient "
                           "contraction here, or ask the density through a gradient-capable face.");
}

// The cached Overlap() accessor is typed for the (symmetric) double cache, which complex bypasses
// (the smat_t->hmat_t cache refactor is pinned).  Use the uncached virtual MakeOverlap() directly --
// for plane waves it is the identity, so this is just trace(D) = the electron count.
template <> double IrrepCD<dcmplx>::GetTotalCharge() const
{
    return std::real(blazem::sum(itsDensityMatrix % blazem::trans(itsBasisSet->MakeOverlap())));
}

template class IrrepCD_Fourier<IrrepCD<dcmplx>>;
template class IrrepCD<dcmplx>;




rvec3_t GradientContraction(const vec_t<rvec3_t>& g, const vec_t<double>& v, const rsmat_t& m)
{
    // No UT coverage
    assert(v.size      ()==m.columns());
    assert(g.size      ()==m.columns());

    rvec3_t ret(0,0,0);
    for (unsigned int i=0; i<v.size(); i++)
        for (unsigned int j=0; j<v.size(); j++)
            ret+=m(i,j)*(g[i]*v[j]+v[i]*g[j]);
    return ret;
}

// rvec3_t GradientContraction(const Vector<Vec3>& g, const Vector<std::complex<double> >& v, const smat_t<std::complex<double> >& m)
// {
//     // No UT coverage
//     assert(v.size      ()==m.columns());
//     assert(g.size      ()==m.columns());

//     Vec3 ret(0,0,0);
//     for (unsigned int i=1; i<=v.size(); i++)
//         for (unsigned int j=1; j<=v.size(); j++)
//             ret+=m(i-1,j-1)*(g(i)*conj(v(j))+v(i)*conj(g(j)));
//     return real(ret);
// }

} //namespace