// File: ExactIrrepCD.C  Exact implementation of the charged density.
module;
#include <cassert>
#include <complex>
#include <iostream>
#include <stdlib.h>
#include <type_traits>
#include <vector>

module qchem.ChargeDensity.Imp.IrrepCD;
import qchem.Symmetry;
import qchem.Blaze;
import qchem.BasisSet.Band_FT_IBS;   // cast the basis UP to the G-space capability (dcmplx path)

namespace qchem::ChargeDensity
{

typedef Vector3D<std::complex<double> > Vec3;

rvec3_t  GradientContraction(const vec_t<rvec3_t >&, const vec_t<double>&, const rsmat_t&);
// rvec3_t  GradientContraction(const Vector<Vec3 >&, const Vector<std::complex<double> >&, const smat_t<std::complex<double> >&);

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> IrrepCD<T>::IrrepCD()
{};

template <class T> IrrepCD<T>::IrrepCD(const DenSMat& D,const tobs_t<T>* theBasisSet,Irrep qns)
    : itsDensityMatrix(D)
    , itsBasisSet(theBasisSet)
    , itsSpin(qns.ms)
    , itsIrrep(qns)
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
template <> void IrrepCD<double>::AccumulateDirect(rsmat_t& Sab, const ohfbs_t* bs_ab) const
{
    const ohfbs_t* bs_cd=dynamic_cast<const ohfbs_t*>(itsBasisSet);
    assert(bs_cd);
    if (!IsZero()) bs_ab->AccumulateDirect(Sab,itsDensityMatrix,bs_cd);
}

template <> void IrrepCD<double>::AccumulateExchange(rsmat_t& Sab, const ohfbs_t* bs_ab) const
{
    const ohfbs_t* bs_cd=dynamic_cast<const ohfbs_t*>(itsBasisSet);
    assert(bs_cd);
    if (!IsZero()) bs_ab->AccumulateExchange(Sab,itsDensityMatrix,bs_cd);
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
// AO density-fit projection <rho|c> = Sum_ab D_ab <ab|c>, the finite (double) path's ProjectedDensity_AO
// face.  The periodic (dcmplx) density is NOT a ProjectedDensity_AO (see ProjectedDensityBase), so this is
// never reached for dcmplx; the if-constexpr keeps the double-only 3-centre machinery out of that build.
template <class T> rvec_t IrrepCD<T>::GetRepulsion3C(const BasisSet::FIT_CD_ABS* fbs) const
{
    if constexpr (std::is_same_v<T,double>)
    {
        if (IsZero()) return rvec_t(fbs->GetNumFunctions(),0.0);
        auto dftbs=dynamic_cast<const todftbs_t<T>*>(itsBasisSet);
        assert(dftbs);
        return dftbs->Repulsion3C(itsDensityMatrix,fbs);
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
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(v->GetMatrix(itsBasisSet,itsSpin,cd)));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

template <class T> double IrrepCD<T>::GetTotalCharge() const
{
    return std::real(blazem::sum(itsDensityMatrix%itsBasisSet->Overlap())); //% is the blaze op for the Shur (direct) product.
}


//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void IrrepCD<T>::ReScale(double factor)
{
    // No UT coverage
    itsDensityMatrix*=factor;
}

template <class T> void IrrepCD<T>::MixIn(const tDM_CD<T>& cd,double c)
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    itsDensityMatrix = itsDensityMatrix*(1-c) + eicd->itsDensityMatrix*c;
}

template <class T> double IrrepCD<T>::GetChangeFrom(const tDM_CD<T>& cd) const
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
template <class T> FourierMap IrrepCD<T>::GetFourierDensity() const
{
    if constexpr (std::is_same_v<T,dcmplx>)
    {
        auto* fb=dynamic_cast<const BasisSet::Band_FT_IBS*>(itsBasisSet);
        assert(fb && "GetFourierDensity requires a Band_FT_IBS (plane-wave) basis");
        return fb->MakeFourierDensity(itsDensityMatrix);
    }
    else
    {
        assert(false && "a finite (non-periodic) density has no reciprocal-lattice Fourier series");
        return FourierMap{};
    }
}


//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& IrrepCD<T>::Write(std::ostream& os) const
{
    return os << itsDensityMatrix;
}

template class IrrepCD<double>;

// --- Complex (plane-wave) density.  HF does NOT apply (the plane-wave path uses no 4-centre exchange), so
// those are NA; the gradient (GGA/plotting) is not yet wired for complex and is unused by the LDA SCF.  The
// AO density-fit projection (GetRepulsion3C) is no longer forced on this path -- IrrepCD<dcmplx> is not a
// ProjectedDensity_AO -- so its old NA-assert specialization is gone.  The remaining members are the generic
// templates above (DM_Contract, op()(r), GetTotalCharge, MixIn, GetChangeFrom, ...), complex-correct.
template <> void IrrepCD<dcmplx>::AccumulateDirect(hmat_t<dcmplx>&, const ohfbs_t*) const
{ assert(false && "AccumulateDirect: HF not applicable to a complex plane-wave density"); }
template <> void IrrepCD<dcmplx>::AccumulateExchange(hmat_t<dcmplx>&, const ohfbs_t*) const
{ assert(false && "AccumulateExchange: HF not applicable to a complex plane-wave density"); }
template <> rvec3_t IrrepCD<dcmplx>::Gradient(const rvec3_t&) const
{ return rvec3_t(0,0,0); }   // No UT coverage; GGA/plotting gradient not yet wired for complex.

// The cached Overlap() accessor is typed for the (symmetric) double cache, which complex bypasses
// (the smat_t->hmat_t cache refactor is pinned).  Use the uncached virtual MakeOverlap() directly --
// for plane waves it is the identity, so this is just trace(D) = the electron count.
template <> double IrrepCD<dcmplx>::GetTotalCharge() const
{
    return std::real(blazem::sum(itsDensityMatrix % blazem::trans(itsBasisSet->MakeOverlap())));
}

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