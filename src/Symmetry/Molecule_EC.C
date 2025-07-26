// File: Symmetry/Molecule_EC.C  Electron configuration for a Molecule.
export module qchem.Symmetry.MoleculeEC;
export import qchem.Symmetry.ElectronConfiguration;

export class Molecule_EC : public virtual ElectronConfiguration
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
