// File: ElectronConfiguration.H
#ifndef _ElectronConfiguration_H_
#define _ElectronConfiguration_H_

class Spin;
class Symmetry;

class ElectronConfiguration
{
public:
    //! Total number of electrons    
    virtual int GetN(                                 ) const=0;
    //! Number of electrons in the spin up/dn sector 
    virtual int GetN(                      const Spin&) const=0;
    //! Number of electrons in one symmetry/irrep sector 
    virtual int GetN(const Symmetry&             ) const=0;
    //! Number of electrons in one {symmetry/irrep,spin up/dn} sector
    virtual int GetN(const Symmetry&, const Spin&) const=0;

    virtual void Display() const=0;
};


const int Nshell=8;
const int LMax=3;

class AtomElectronConfiguration : public virtual ElectronConfiguration
{
public: 
    AtomElectronConfiguration(int Z);
    //AtomElectronConfiguration(int Z, double charge, int val_config[4]);
    
    virtual int GetN(                                 ) const;
    virtual int GetN(                      const Spin&) const;
    virtual int GetN(const Symmetry&             ) const;
    virtual int GetN(const Symmetry&, const Spin&) const;
    virtual void Display() const;

private:
    static const int FullShells[Nshell][LMax+2];
    int N [LMax+1]; //Total N[l]
    int Nf[LMax+1]; //Full shell N[l]
    int Nv[LMax+1]; //Valance shell N[l]
    int NUnpaired;
    double charge;
};

class MoleculeElectronConfiguration : public virtual ElectronConfiguration
{
public: 
    MoleculeElectronConfiguration() : Ne(0) {};
    MoleculeElectronConfiguration(int _Ne) : Ne(_Ne) {};
    //AtomElectronConfiguration(int Z, double charge, int val_config[4]);
    
    virtual int GetN(                                 ) const {return Ne;}
    virtual int GetN(                      const Spin&) const;
    virtual int GetN(const Symmetry&             ) const {return Ne;}
    virtual int GetN(const Symmetry&, const Spin& s) const {return GetN(s);}
    virtual void Display() const;
private:
    int Ne;
};

#endif // _ElectronConfiguration_H_
