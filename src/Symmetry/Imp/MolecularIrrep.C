// File: Symmetry/Imp/MolecularIrrep.C
module;
#include <ostream>
module qchem.Symmetry.MolecularIrrep;

namespace qchem {

std::ostream& MolecularIrrep::Write(std::ostream& os) const
{
    return os << itsLabel;
}

} // namespace qchem