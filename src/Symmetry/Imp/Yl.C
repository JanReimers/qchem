// File: Symmetry/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iostream>
#include <cassert>

module qchem.Symmetry.Yl;
import qchem.Common.Strings;


Yl_Sym::Yl_Sym(int theL)
    : itsL(theL)
{};

 size_t Yl_Sym::SequenceIndex() const //Used for op<
 {
    return itsL;
 }


int Yl_Sym::GetDegeneracy() const
{
    return 2*itsL+1;
}




std::ostream& Yl_Sym::Write(std::ostream& os) const
{
    return os << SPDFG[itsL] << " ";
}

