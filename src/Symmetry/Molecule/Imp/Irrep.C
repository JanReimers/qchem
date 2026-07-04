// File: Symmetry/Molecule/Imp/Irrep.C
module;
#include <ostream>
module qchem.Symmetry.Molecule.Irrep;

namespace qchem::Symmetry::Molecule {

std::ostream& Irrep::Write(std::ostream& os) const
{
    return os << itsLabel;
}

} // namespace qchem::Symmetry::Molecule