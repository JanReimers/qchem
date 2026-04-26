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

export struct ml_Breakdown
{
    std::vector<int> ml_paired;     //List of ml values for paired orbitals
    std::vector<int> ml_unpaired;   //List of ml values for unpaired orbitals
    std::vector<int> ml_unoccupied; //List of ml values for empty orbitals
};

export class Atom_EC : public virtual ElectronConfiguration
{
public: 
    typedef std::set<Irrep_QNs::sym_t> syms_t;
    Atom_EC(int Z);
    
    virtual int    GetN(const Irrep_QNs&) const;  //Core + Valance
    virtual size_t GetLMax() const {return itsLMax;}
    virtual void   Display() const;
    
    ml_Breakdown GetBreadown(size_t l) const;
    syms_t GetIrreps() const;
private:
    typedef std::shared_ptr<const Symmetry> sym_t;
    friend class ElectronConfigurationTests;
    int  GetN() const;
    int  GetN(const Spin&) const;
    int  GetN(const sym_t&) const;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    double charge;
    size_t itsLMax,itsLValance;
    std::map<Irrep_QNs,size_t> itsOccupations; //Spin polarized list;
    std::map<Irrep_QNs,size_t> itsUnpolOccupations; //Spin un polarized list;
};

