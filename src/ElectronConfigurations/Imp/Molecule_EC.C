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
    if (qns.ms==Spin::None) return Ne;
    if (Ne%2==0)
        return Ne/2;
    else
        return qns.ms==Spin::Up ? (Ne+1)/2 : (Ne-1)/2;
}

ElectronConfiguration::syms_t Molecule_EC::GetIrreps() const
{
    syms_t ret;
    ret.insert(Symmetry::UnitFactory());
    return ret;
}
void Molecule_EC::Display() const
{
    cout << "Ne: " << Ne << endl;
}

} // namespace qchem