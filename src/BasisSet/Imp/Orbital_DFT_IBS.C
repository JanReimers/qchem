// File: BasisSet/Imp/Orbital_DFT_IBS.C
module;
#include <cassert>
module qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace qchem::BasisSet
{
// The tensor lives in theCache<TFit>() (NOT <T>): its element type follows the FIT axis, so that is the
// cache that can store it.  Identical for every T==TFit instantiation; load-bearing for <double,dcmplx>
// (a real TRIM block sharing the run's complex G-space fit basis -- doc/RealComplexPlan.md).
template <class T, class TFit> const Projector3<TFit>& Orbital_DFT_IBS<T,TFit>::Overlap3C  (const FIT_SF_ABS<TFit>& c) const
{
    return theCache<TFit>().Get(IntegralsCache_Base::I3C::Overlap,this,&c,
        [this,&c]{ return MakeOverlap3C(c); });
}
//! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
template <class T, class TFit> const Projector3<TFit>& Orbital_DFT_IBS<T,TFit>::Repulsion3C(const FIT_CD_ABS<TFit>& c) const
{
    return theCache<TFit>().Get(IntegralsCache_Base::I3C::Repulsion,this,&c,
        [this,&c]{ return MakeRepulsion3C(c); });
}

// The density-matrix contraction <rho|c> = Sum_ab D_ab <ab|c> that used to live here (as the convenience
// overloads Overlap3C(D,c)/Repulsion3C(D,c)) now lives in the charge density (IrrepCD::GetRepulsion3C),
// which owns D -- the basis exposes only the D-free integral tensors Overlap3C(c)/Repulsion3C(c) above.

template class Orbital_DFT_IBS<double>;
template class Orbital_DFT_IBS<dcmplx>;   // the reciprocal-space lineage (ex Band_FT_IBS -- V1.1)
}
