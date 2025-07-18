// File: Molecule_EC.H
#ifndef _Molecule_EC_H_
#define _Molecule_EC_H_

#include <Symmetry/ElectronConfiguration.H>
import qchem.Symmetry.Spin;
import qchem.Symmetry;

class Molecule_EC : public virtual ElectronConfiguration
{
public: 
    Molecule_EC() : Ne(0) {};
    Molecule_EC(int _Ne) : Ne(_Ne) {};
    
    virtual int GetN(const Irrep_QNs& qns) const;
    virtual void Display() const;
private:
    int GetN() const {return Ne;}
    int GetN(const Spin&) const;
    int GetN(const Symmetry&) const {return Ne;}
    int Ne;
};

#endif // _Molecule_EC_H_
