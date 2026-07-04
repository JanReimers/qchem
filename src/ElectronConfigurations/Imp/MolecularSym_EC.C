// File: ElectronConfigurations/Imp/MolecularSym_EC.C
module;
#include <map>
#include <string>
#include <iostream>
module qchem.ElectronConfiguration.MolecularSym;
import qchem.Symmetry.Molecule.Irrep;

namespace qchem {

int MolecularSym_EC::GetN(const Irrep& qns) const
{
    std::string label = qns.sym ? qns.sym->GetLabel() : "";
    auto it = itsOcc.find(label);
    int n = (it!=itsOcc.end()) ? it->second : 0;
    if (qns.ms==Spin::None) return n;            // unpolarized: both spins
    return (qns.ms==Spin::Up) ? (n+1)/2 : n/2;   // polarized split
}

ElectronConfiguration::syms_t MolecularSym_EC::GetIrreps() const
{
    syms_t ret;
    int i=0;
    for (const auto& [label,n] : itsOcc) ret.insert(sym_t(new Symmetry::Molecule::Irrep(label, i++)));
    return ret;
}

void MolecularSym_EC::Display() const
{
    for (const auto& [label,n] : itsOcc) std::cout << label << ":" << n << "  ";
    std::cout << std::endl;
}

} // namespace qchem