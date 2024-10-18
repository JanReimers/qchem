// File: Orbital.C
#include "Orbital.H"

template class TOrbital<double>;
template class TOrbital<std::complex<double> >;
template class TOrbitals<double>;
template class TOrbitals<std::complex<double> >;

//#include "BasisSet.H"
//#include "Misc/rc_ptr.H"
//#include <cassert>

//void OrbitalGroup::FixUpPointer(const rc_ptr<const IrrepBasisSet>& bs)
//{
//    // No UT coverage
//    assert(&*bs);
//    for (auto o:*this) o->FixUpPointer(&*bs);
//}

