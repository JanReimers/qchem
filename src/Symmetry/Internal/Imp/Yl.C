// File: Symmetry/Internal/Imp/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iostream>

module qchem.Symmetry.Internal.Spherical;
import qchem.Common.Strings;

namespace Symmetryns::Internal::Spherical
{

std::ostream& Yl_Sym::Write(std::ostream& os) const
{
    return os << SPDFG[itsL] << " ";
}

} //namespace
