// File: Atom_EC.C  Electron configuration for atoms.
module;
#include <vector>
#include <set>
#include <map>
#include <memory>
#include "forward.H"

export module qchem.Symmetry.AtomEC;
import qchem.Symmetry;
import qchem.Symmetry.ElectronCounts;
import qchem.Symmetry.ElectronConfiguration;
const int Nshell=8;
using namespace Symmetry;

export class Atom_EC : public virtual ElectronConfiguration
{
public: 
    static constexpr auto cmp = [](sym_t a, sym_t b) 
        {
             return a->SequenceIndex() < b->SequenceIndex();
        }; 
    typedef std::set<sym_t,decltype(cmp)> syms_t;
    Atom_EC(int Z);
    
    virtual int    GetN(const Irrep_QNs&) const;  //Core + Valance
    virtual size_t GetLMax() const {return itsLMax;}
    virtual void   Display() const;
    
    syms_t GetIrreps() const;
protected:
    friend class ElectronConfigurationTests;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    size_t itsLMax,itsLValance;
    std::map<Irrep_QNs,size_t> itsOccupations; //Spin polarized list;
    std::map<Irrep_QNs,size_t> itsUnpolOccupations; //Spin un polarized list;
};

