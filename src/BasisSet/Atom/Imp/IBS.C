// File: BasisSet/Atom/Imp/IBS.C Atom specific irrep basis sets.
module;
#include <iostream>
module qchem.BasisSet.Atom.IBS;
import qchem.Streamable;

namespace AtomBS
{
std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << itsEval->Name() << " symmetry=" << GetSymmetry();
    itsEval->Write(os);
    return os;
}

}//namespace
