// File: Atom_EC.H
#ifndef _Atom_EC_H_
#define _Atom_EC_H_

#include <ElectronConfiguration.H>

const int Nshell=8;
const int LMax=3;

struct ElCounts
{
    ElCounts() : N{0,0,0,0} {};
    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell N[l]
    int Nv[LMax+1]; //Valance shell N[l]
    int NUnpaired;
};

class Atom_EC : public virtual ElectronConfiguration
{
public: 
    Atom_EC(int Z);
    
    virtual int    GetN(const Irrep_QNs&) const;
    virtual size_t GetLMax() const {return itsLMax;}
    virtual void   Display() const;

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
