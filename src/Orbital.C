// File: Orbital.C
#include "Orbital.H"

template class TOrbital<double>;
template class TOrbital<std::complex<double> >;
template class TOrbitalGroup<double>;
template class TOrbitalGroup<std::complex<double> >;

#include "BasisSet.H"
#include "Misc/rc_ptr.H"
#include <cassert>

void OrbitalGroup::FixUpPointer(const rc_ptr<const BasisSet>& bs)
{
    // No UT coverage
    assert(&*bs);
    for (auto o:*this) o->FixUpPointer(&*bs);
}

