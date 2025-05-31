// File: Atom_EC.H
#ifndef _Atom_EC_H_
#define _Atom_EC_H_

#include <ElectronConfiguration.H>
#include <vector>

const int Nshell=8;
const int LMax=3;

struct ElCounts
{
    ElCounts() : N{0,0,0,0}, Nf{0,0,0,0},Nv{0,0,0,0},Nu{0,0,0,0} {};
    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell N[l]
    int Nv[LMax+1]; //Valance shell N[l]
    int Nu[LMax+1]; //# of unpaired electrons for a given l.
};

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
    
    virtual int    GetN(const Irrep_QNs&) const;
    virtual size_t GetLMax() const {return itsLMax;}
    virtual size_t GetNUnapired(size_t l) const {return itsNs.Nu[l];}
    virtual size_t GetNval(size_t l) const {return itsNs.Nv[l];}
    virtual void   Display() const;
    
    ml_Breakdown GetBreadown(size_t l) const;

private:
    friend class ElectronConfigurationTests;
    int  GetN() const;
    int  GetN(const Spin&) const;
    int  GetN(const Symmetry&) const;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    double charge;
    size_t itsLMax;
};

#endif // _Atom_EC_H_
