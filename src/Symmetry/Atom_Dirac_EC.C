// File: Symmetry/Atom_Dirac_EC.C  Electron configuration for Relatistic/Dirac atoms.
module;
#include "forward.H"

export module qchem.Symmetry.Atom_Dirac_EC;
import qchem.Symmetry.AtomEC;

export class Atom_Dirac_EC 
    : public virtual ElectronConfiguration
    , private Atom_EC
{
public: 
    using Atom_EC::GetN;
    using Atom_EC::GetIrreps;
    Atom_Dirac_EC(int Z);
    
private:
    friend class ElectronConfigurationTests;
};

