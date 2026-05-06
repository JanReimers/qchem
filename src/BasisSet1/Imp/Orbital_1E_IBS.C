// File: BasisSet1/Imp/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
#include <cassert>
#include <string>
module qchem.BasisSet1.Orbital_1E_IBS;

import qchem.BasisSet1.DB_Cache;

namespace BasisSet1
{

template <class T> const smat_t<T>& Integrals_Kinetic<T>::Kinetic() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Has(IntegralsCache_Base::I2C::Kinetic,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
        ? cache->GetSMat() : cache->Set(MakeKinetic());
    
}

template <class T> const smat_t<T>& Integrals_Nuclear<T>::Nuclear(const Cluster* cl) const
{
    assert(cl);
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Has(IntegralsCache_Base::I2n::Nuclear,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),cl->ID())
        ? cache->GetSMat() : cache->Set(MakeNuclear(cl));
}

template class Integrals_Kinetic<double>;
template class Integrals_Nuclear<double>;

} //namespace 