// File: ElectronConfigurations/Molecule_EC.C  Electron configuration for a Molecule.
export module qchem.ElectronConfiguration.Molecule;
export import qchem.ElectronConfiguration;

export class Molecule_EC : public virtual ElectronConfiguration
{
public: 
    Molecule_EC() : Ne(0) {};
    Molecule_EC(int _Ne) : Ne(_Ne) {};
    
    virtual int    GetN(const Irrep& qns) const;
    virtual syms_t GetIrreps() const;
    virtual void Display() const;
private:
    int Ne;
};
