// File: BasisSet/Imp/Orbital_DFT_IBS.C
module;
#include <memory>
#include <cassert>
module qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.BasisCache;

namespace qchem::BasisSet
{
// The tensor lives in the PARENT's cache, keyed by (this block's index, the fit's INSTANCE serial) (doc/IntegralCacheOwnershipPlan.md
// S3a, §3b): never the fit's address (reused after it dies) nor its content ID.  Its element type follows the FIT axis (TFit, NOT T),
// so it lands in the cache's TFit map.  Identical for every T==TFit instantiation; load-bearing for <double,dcmplx>
// (a real TRIM block sharing the run's complex G-space fit basis -- doc/RealComplexPlan.md).
template <class T, class TFit> const PairProjection<TFit>& Orbital_DFT_IBS<T,TFit>::Overlap3C  (const FIT_SF_ABS<TFit>& c) const
{
    return this->Cache().template Get<TFit>(BasisCache::Three::Overlap, this->CacheIndex(), c.Instance(), c.Name(),
        [this,&c]{ return MakeOverlap3C(c); });
}
//! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
template <class T, class TFit> const PairProjection<TFit>& Orbital_DFT_IBS<T,TFit>::Repulsion3C(const FIT_CD_ABS<TFit>& c) const
{
    return this->Cache().template Get<TFit>(BasisCache::Three::Repulsion, this->CacheIndex(), c.Instance(), c.Name(),
        [this,&c]{ return MakeRepulsion3C(c); });
}

// The density-matrix contraction <rho|c> = Sum_ab D_ab <ab|c> that used to live here (as the convenience
// overloads Overlap3C(D,c)/Repulsion3C(D,c)) now lives in the charge density (IrrepCD::GetRepulsion3C),
// which owns D -- the basis exposes only the D-free integral tensors Overlap3C(c)/Repulsion3C(c) above.

template class Orbital_DFT_IBS<double>;
template class Orbital_DFT_IBS<dcmplx>;         // the reciprocal-space lineage (ex Band_FT_IBS -- V1.1)
template class Orbital_DFT_IBS<double,dcmplx>;  // a REAL TRIM block on the run's complex fit basis (Step 3)
}
