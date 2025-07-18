// File: ElectronConfiguration.C
export module qchem.Symmetry.ElectronConfiguration;
export import qchem.Symmetry.Irrep;

export class ElectronConfiguration
{
public:
    virtual int    GetN(const Irrep_QNs&) const=0;
    virtual void   Display() const=0;
};