// File: ElectronConfigurations/MolecularSym_EC.C
// Electron configuration for a symmetric molecule with a KNOWN per-irrep occupation, keyed by
// Mulliken irrep label (e.g. water C2v: {A1:6, B1:2, B2:2}).  A stop-gap for the first
// symmetric-HF integration test; the general case (global aufbau across irreps) belongs in a
// CompositeAufbauWF that overrides FillOrbitals.
module;
#include <map>
#include <string>
export module qchem.ElectronConfiguration.MolecularSym;
export import qchem.ElectronConfiguration;

export class MolecularSym_EC : public virtual ElectronConfiguration
{
public:
    // occ: total electrons per Mulliken irrep label, e.g. {{"A1",6},{"B1",2},{"B2",2}}.
    MolecularSym_EC(const std::map<std::string,int>& occ) : itsOcc(occ) {}

    virtual int    GetN(const Irrep& qns) const;
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;

private:
    std::map<std::string,int> itsOcc;
};
