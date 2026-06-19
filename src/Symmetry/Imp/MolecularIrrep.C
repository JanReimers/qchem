// File: Symmetry/Imp/MolecularIrrep.C
module;
#include <ostream>
module qchem.Symmetry.MolecularIrrep;

std::ostream& MolecularIrrep::Write(std::ostream& os) const
{
    return os << itsLabel;
}
