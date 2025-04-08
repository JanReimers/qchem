// File: ElectronConfiguration.H
#ifndef _ElectronConfiguration_H_
#define _ElectronConfiguration_H_

class Spin;
class Symmetry;

class ElectronConfiguration
{
public:
    //! Total number of electrons    
    virtual int GetN(                            ) const=0;
    //! Number of electrons in the spin up/dn sector 
    virtual int GetN(                 const Spin&) const=0;
    //! Number of electrons in one symmetry/irrep sector 
    virtual int GetN(const Symmetry&             ) const=0;
    //! Number of electrons in one {symmetry/irrep,spin up/dn} sector
    virtual int GetN(const Symmetry&, const Spin&) const=0;

    virtual void Display() const=0;
};


const int Nshell=8;
const int LMax=3;

class Atom_EC : public virtual ElectronConfiguration
{
public: 
    Atom_EC(int Z);
    
    virtual int  GetN(                            ) const;
    virtual int  GetN(                 const Spin&) const;
    virtual int  GetN(const Symmetry&             ) const;
    virtual int  GetN(const Symmetry&, const Spin&) const;
    virtual void Display() const;

private:
    static const int FullShells[Nshell][LMax+2];
    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell N[l]
    int Nv[LMax+1]; //Valance shell N[l]
    int NUnpaired;
    double charge;
};

class Molecule_EC : public virtual ElectronConfiguration
{
public: 
    Molecule_EC() : Ne(0) {};
    Molecule_EC(int _Ne) : Ne(_Ne) {};
    
    virtual int GetN(                              ) const {return Ne;}
    virtual int GetN(                   const Spin&) const;
    virtual int GetN(const Symmetry&               ) const {return Ne;}
    virtual int GetN(const Symmetry&, const Spin& s) const {return GetN(s);}
    virtual void Display() const;
private:
    int Ne;
};

#endif // _ElectronConfiguration_H_
