// File: ElectronConfiguration.H
#ifndef _ElectronConfiguration_H_
#define _ElectronConfiguration_H_

#include <types.H>

class ElectronConfiguration
{
public:
    virtual int GetN(const Irrep_QNs&) const=0;
    virtual void Display() const=0;
};


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
    
    virtual int  GetN(const Irrep_QNs&) const;
    virtual void Display() const;

private:
    friend class ElectronConfigurationTests;
    int  GetN() const;
    int  GetN(const Spin&) const;
    int  GetN(const Symmetry&) const;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    double charge;
};

class Molecule_EC : public virtual ElectronConfiguration
{
public: 
    Molecule_EC() : Ne(0) {};
    Molecule_EC(int _Ne) : Ne(_Ne) {};
    
    virtual int GetN(const Irrep_QNs& qns) const;
    virtual void Display() const;
private:
    int GetN(                              ) const {return Ne;}
    int GetN(                   const Spin&) const;
    int GetN(const Symmetry&               ) const {return Ne;}
    int Ne;
};

#endif // _ElectronConfiguration_H_
