// File: Symmetry/ElectronConfiguration.C Interface for and electron configuration.
export module qchem.Symmetry.ElectronConfiguration;
export import qchem.Symmetry.Irrep;

export class ElectronConfiguration
{
public:
    virtual int    GetN(const Irrep&) const=0;
    virtual void   Display() const=0;
};