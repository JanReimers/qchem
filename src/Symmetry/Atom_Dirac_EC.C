// File: Symmetry/Atom_Dirac_EC.C  Electron configuration for Relatistic/Dirac atoms.
module;
#include <vector>
#include <set>
#include <map>
#include <memory>
#include "forward.H"

export module qchem.Symmetry.Atom_Dirac_EC;
import qchem.Symmetry.AtomEC;
import qchem.Symmetry;
import qchem.Symmetry.ElectronCounts;
import qchem.Symmetry.ElectronConfiguration;
const int Nshell=8;

export class Atom_Dirac_EC 
    : public virtual ElectronConfiguration
    , private Atom_EC
{
public: 
    using Atom_EC::cmp;
    using Atom_EC::syms_t;
    using Atom_EC::GetN;
    using Atom_EC::GetIrreps;
    Atom_Dirac_EC(int Z) : Atom_EC(Z) {};
    
    // virtual int    GetN(const Irrep_QNs&) const;  //Core + Valance
    // virtual size_t GetLMax() const {return itsLMax;}
    // virtual void   Display() const;
    
    // syms_t GetIrreps() const;
private:
    friend class ElectronConfigurationTests;
    // int  GetN() const;
    // int  GetN(const Spin&) const;
    // int  GetN(const Irrep_QNs::sym_t&) const;

    // static const int FullShells[Nshell][LMax+2];
    // ElCounts itsNs; //Total,core, valance and unpaired counts.
    // size_t itsLMax,itsLValance;
    // std::map<Irrep_QNs,size_t> itsOccupations; //Spin polarized list;
    // std::map<Irrep_QNs,size_t> itsUnpolOccupations; //Spin un polarized list;
};

