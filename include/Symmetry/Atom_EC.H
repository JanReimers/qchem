// File: Atom_EC.H
#ifndef _Atom_EC_H_
#define _Atom_EC_H_

#include "Symmetry/ElectronCounts.H"
#include <Symmetry/ElectronConfiguration.H>
#include <vector>
#include <memory>

const int Nshell=8;

struct ml_Breakdown
{
    std::vector<int> ml_paired;     //List of ml values for paired orbitals
    std::vector<int> ml_unpaired;   //List of ml values for unpaired orbitals
    std::vector<int> ml_unoccupied; //List of ml values for empty orbitals
};

class Atom_EC : public virtual ElectronConfiguration
{
public: 
    Atom_EC(int Z);
    
    virtual int    GetN(const Irrep_QNs&) const;  //Core + Valance
    virtual size_t GetLMax() const {return itsLMax;}
    virtual void   Display() const;
    virtual bool   IsMagnetic() const; //Does this EC require ml splitting for the basis set?
    
    ml_Breakdown GetBreadown(size_t l) const;

private:
    typedef std::shared_ptr<const Symmetry> sym_t;
    friend class ElectronConfigurationTests;
    int  GetN() const;
    int  GetN(const Spin&) const;
    int  GetN(const sym_t&) const;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    double charge;
    size_t itsLMax;
};

#endif // _Atom_EC_H_
