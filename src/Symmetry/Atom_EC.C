// File: Symmetry/Atom_EC.C  Electron configuration for atoms.
module;
#include <map>
#include "forward.H"

export module qchem.Symmetry.AtomEC;
export import qchem.Symmetry.ElectronConfiguration;

import qchem.Symmetry.ElectronCounts;
const int Nshell=8;
// using namespace Symmetry;

export class Atom_EC : public virtual ElectronConfiguration
{
public:  
    Atom_EC(int Z);
    virtual int    GetN(const Irrep&) const;  //Core + Valance
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;  

    size_t GetLMax() const {return itsLMax;}
protected:
    friend class ElectronConfigurationTests;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    size_t itsLMax,itsLValance;
    std::map<Irrep,size_t> itsOccupations; //Spin polarized list;
    std::map<Irrep,size_t> itsUnpolOccupations; //Spin un polarized list;
};

