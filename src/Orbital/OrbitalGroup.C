// File: OrbitalGroup.C  Group of orbitals.



#include "Orbital/Orbital.H"
#include "Orbital/OrbitalGroup.H"
//#include "Orbital/OrbitalGroupBrowser.H"
#include "BasisSet/BasisSet.H"
#include "Misc/rc_ptr.H"
#include <cassert>

void OrbitalGroup::FixUpPointer(const rc_ptr<const BasisSet>& bs)
{
    // No UT coverage
    assert(&*bs);
//    for (OrbitalGroupIterator i(*this); i; i++) i->FixUpPointer(&*bs);
    for (auto o:*this) o->FixUpPointer(&*bs);
}


