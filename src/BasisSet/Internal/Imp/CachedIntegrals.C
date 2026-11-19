// File: BasisSet/Internal/Imp/CachedIntegrals.C  The cached kinetic and 3-centre accessors.
module;
module qchem.BasisSet.Internal.CachedIntegrals;

namespace qchem::BasisSet
{

template <class T> const hmat_t<T>& CachedKinetic<T>::Kinetic() const
{
    return this->Cache().template Get<T>(BasisCache::Own::Kinetic, this->CacheIndex(), [this]{ return this->MakeKinetic(); });
}
template class CachedKinetic<double>;
template class CachedKinetic<dcmplx>;   // the periodic 1E bases (plane waves, GPW)

template <class T, class TFit> const PairProjection<TFit>& CachedFit3C<T,TFit>::Overlap3C(const FIT_SF_ABS<TFit>& c) const
{
    return this->Cache().template Get<TFit>(BasisCache::Three::Overlap, this->CacheIndex(), c.Instance(), c.Name(),
        [this,&c]{ return this->MakeOverlap3C(c); });
}
//! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
template <class T, class TFit> const PairProjection<TFit>& CachedFit3C<T,TFit>::Repulsion3C(const FIT_CD_ABS<TFit>& c) const
{
    return this->Cache().template Get<TFit>(BasisCache::Three::Repulsion, this->CacheIndex(), c.Instance(), c.Name(),
        [this,&c]{ return this->MakeRepulsion3C(c); });
}
template class CachedFit3C<double>;
template class CachedFit3C<dcmplx>;          // the reciprocal-space lineage
template class CachedFit3C<double,dcmplx>;   // a REAL TRIM block on the run's complex fit basis (doc/RealComplexPlan.md)

} //namespace
