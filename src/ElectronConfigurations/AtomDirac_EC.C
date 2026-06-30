// File: ElectronConfigurations/AtomDirac_EC.C  Electron configuration for Relatistic/Dirac atoms.
module;
#include "forward.H"

export module qchem.ElectronConfiguration.AtomDirac;
import qchem.ElectronConfiguration.AtomNR;

namespace qchem {

export class AtomDirac_EC 
    : public virtual ElectronConfiguration
    , private Atom_EC
{
public: 
    using Atom_EC::GetN;
    using Atom_EC::GetIrreps;
    AtomDirac_EC(int Z);
    
private:
    friend class ::ElectronConfigurationTests;
};


} // namespace qchem