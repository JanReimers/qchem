// File: OrbitalGroup.C  Group of orbitals.



#include "Orbital/Orbital.H"
#include "Orbital/OrbitalGroup.H"
#include "BasisSet/BasisSet.H"
#include "Misc/rc_ptr.H"
#include <cassert>

void OrbitalGroup::FixUpPointer(const rc_ptr<const BasisSet>& bs)
{
    // No UT coverage
    assert(&*bs);
    for (auto o:*this) o->FixUpPointer(&*bs);
}


