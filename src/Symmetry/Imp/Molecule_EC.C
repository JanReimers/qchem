// File: Symmetry/Molecule_ECImp.C  Electron configuration for a Molecule.
module;
#include <cassert>
#include <iostream>
module qchem.Symmetry.MoleculeEC;
import qchem.Symmetry.Irrep;

using std::cout;
using std::endl;

int Molecule_EC::GetN(const Irrep_QNs& qns) const
{
    return GetN(qns.ms);
}
int Molecule_EC::GetN(const Spin& s) const
{
    if (s==Spin::None) return GetN();
    if (Ne%2==0)
        return Ne/2;
    else
        return s==Spin::Up ? (Ne+1)/2 : (Ne-1)/2;
}

void Molecule_EC::Display() const
{
    cout << "Ne: " << Ne << endl;
}
