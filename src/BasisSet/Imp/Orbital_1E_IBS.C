// File: BasisSet/Imp/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
#include <memory>
#include <cassert>
module qchem.BasisSet.Orbital_1E_IBS;

import qchem.BasisSet.BasisCache;

namespace qchem::BasisSet
{

// Cached accessor for the kinetic BUILDING BLOCK \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$
// (NOT the kinetic energy -- no 1/2; see BasisSet/Orbital_1E_IBS.C).
template <class T> const hmat_t<T>& Integrals_Kinetic<T>::Kinetic() const
{
    return this->Cache().template Get<T>(BasisCache::Own::Kinetic, this->CacheIndex(), [this]{ return MakeKinetic(); });
}

template class Integrals_Kinetic<double>;
// The complex flavours: the periodic (dcmplx) 1E bases -- plane waves and now GPW -- use these cached
// accessors too (mirrors Integrals_Overlap<dcmplx> in Imp/IrrepBasisSet.C).
template class Integrals_Kinetic<dcmplx>;

} //namespace