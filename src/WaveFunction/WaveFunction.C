// File: WaveFunction.C  Interface for a wave function.



#include "WaveFunction/WaveFunction.H"
#include "Orbital/TOrbitalGroup.H"
#include <cassert>

void WaveFunction::FixUpPointer(OrbitalGroup* og, const rc_ptr<const BasisSet>& bs)
{
    assert(og);
    og->FixUpPointer(bs);
}

