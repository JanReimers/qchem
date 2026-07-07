// File: BasisSet/Imp/Orbital_DFT_IBS.C
module;
#include <cassert>
module qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace qchem::BasisSet
{
template <class T> const ERI3<T>& Orbital_DFT_IBS<T>::Overlap3C  (const rFIT_SF_ABS& c) const
{
    return theCache<T>().Get(IntegralsCache_Base::I3C::Overlap,this,&c,
        [this,&c]{ return MakeOverlap3C(c); });
} 
//! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
template <class T> const ERI3<T>& Orbital_DFT_IBS<T>::Repulsion3C(const rFIT_CD_ABS& c) const
{
    return theCache<T>().Get(IntegralsCache_Base::I3C::Repulsion,this,&c,
        [this,&c]{ return MakeRepulsion3C(c); });
} 

// The density-matrix contraction <rho|c> = Sum_ab D_ab <ab|c> that used to live here (as the convenience
// overloads Overlap3C(D,c)/Repulsion3C(D,c)) now lives in the charge density (IrrepCD::GetRepulsion3C),
// which owns D -- the basis exposes only the D-free integral tensors Overlap3C(c)/Repulsion3C(c) above.

template class Orbital_DFT_IBS<double>;
}
