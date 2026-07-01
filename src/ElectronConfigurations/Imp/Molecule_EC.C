// File: ElectronConfigurations/Imp/Molecule_EC.C  Electron configuration for a Molecule.
module;
#include <iostream>
module qchem.ElectronConfiguration.Molecule;
import qchem.Symmetry.Factory;

namespace qchem {
using std::cout;
using std::endl;

int Molecule_EC::GetN(const Irrep& qns) const
{
    if (qns.ms==Spin::Up)   return itsNup;
    if (qns.ms==Spin::Down) return itsNdn;
    return itsNup + itsNdn;   // Spin::None -> total (unpolarized) count
}

ElectronConfiguration::syms_t Molecule_EC::GetIrreps() const
{
    syms_t ret;
    ret.insert(Symmetry::UnitFactory());
    return ret;
}
void Molecule_EC::Display() const
{
    cout << "Nup: " << itsNup << " Ndn: " << itsNdn << " (Ne=" << itsNup+itsNdn << ")" << endl;
}

} // namespace qchem