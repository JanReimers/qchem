// File: ElectronConfigurations/Molecule_EC.C  Electron configuration for a Molecule.
export module qchem.ElectronConfiguration.Molecule;
export import qchem.ElectronConfiguration;

namespace qchem {

//! Molecular electron configuration.  SPIN-NATIVE: the primary form is the two channel counts
//! (nUp, nDown); the single-argument Molecule_EC(Ne) is the minimal-spin collapse (singlet for even Ne,
//! doublet for odd) -- the closed-shell convenience, NOT the default formulation.  Aufbau fills each spin
//! channel independently across the point-group irreps (UsesAufbau), so open-shell (nUp != nDown) needs no
//! SCF-loop change -- it just reports asymmetric per-spin totals.  Front-end multiplicity 2S+1 maps here as
//! nUp-nDown = 2S, nUp+nDown = Ne (see the facade, OpenWork B4).
export class Molecule_EC : public virtual ElectronConfiguration
{
public:
    Molecule_EC() : itsNup(0), itsNdn(0) {};
    //! Closed-shell / minimal-spin: Ne electrons, smallest spin (the zeta=0 collapse of the (nUp,nDown) form).
    Molecule_EC(int Ne) : itsNup((Ne+1)/2), itsNdn(Ne/2) {};
    //! Spin-native open-shell: \a nUp spin-up and \a nDown spin-down electrons (the magnetism path).
    Molecule_EC(int nUp, int nDown) : itsNup(nUp), itsNdn(nDown) {};

    virtual int    GetN(const Irrep& qns) const;
    virtual syms_t GetIrreps() const;
    virtual void Display() const;
    virtual bool UsesAufbau() const {return true;}   // molecular aufbau across point-group irreps
private:
    int itsNup, itsNdn;
};

} // namespace qchem