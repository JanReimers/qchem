// File: OrbitalGroup.C  Group of orbitals.



#include "Orbital.H"
#include "OrbitalGroup.H"
#include "BasisSet.H"
#include "Misc/rc_ptr.H"
#include <cassert>

void OrbitalGroup::FixUpPointer(const rc_ptr<const BasisSet>& bs)
{
    // No UT coverage
    assert(&*bs);
    for (auto o:*this) o->FixUpPointer(&*bs);
}


