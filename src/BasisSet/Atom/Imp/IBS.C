// File: BasisSet/Atom/Imp/IBS.C Atom specific irrep basis sets.
module;
#include <iostream>
module qchem.BasisSet.Atom.IBS;

namespace Atom
{
std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << "Spherical Slater L=" << *GetSymmetry();
    itsEval->Write(os);
    return os;
}

}//namespace
