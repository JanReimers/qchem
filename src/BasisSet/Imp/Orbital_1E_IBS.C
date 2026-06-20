// File: BasisSet/Imp/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
#include <cassert>
module qchem.BasisSet.Orbital_1E_IBS;

import qchem.BasisSet.Internal.DB_Cache;

namespace BasisSet
{

// Cached accessor for the kinetic BUILDING BLOCK \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$
// (NOT the kinetic energy -- no 1/2; see BasisSet/Orbital_1E_IBS.C).
template <class T> const smat_t<T>& Integrals_Kinetic<T>::Kinetic() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::Kinetic,this,
        [this]{ return MakeKinetic(); });
}

template <class T> const smat_t<T>& Integrals_Nuclear<T>::Nuclear(const Cluster* cl) const
{
    assert(cl);
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2n::Nuclear,this,cl->ID(),
        [this,cl]{ return MakeNuclear(cl); });
}

template class Integrals_Kinetic<double>;
template class Integrals_Nuclear<double>;

} //namespace 