// File: Symmetry/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iostream>
#include <cassert>

module qchem.Symmetry.Yl;
import qchem.Common.Strings;


std::ostream& Yl_Sym::Write(std::ostream& os) const
{
    return os << SPDFG[itsL] << " ";
}

