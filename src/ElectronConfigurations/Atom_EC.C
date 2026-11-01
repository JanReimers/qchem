// File: ElectronConfigurations/Atom_EC.C  Electron configuration for atoms.
module;
#include <map>
#include "forward.H"

export module qchem.ElectronConfiguration.AtomNR;
export import qchem.ElectronConfiguration;

import qchem.ElectronConfiguration.ElectronCounts;

namespace qchem {
const int Nshell=8;

export class Atom_EC : public virtual ElectronConfiguration
{
public:
    Atom_EC(int Z);
    virtual int    GetN(const Irrep&) const;  //Core + Valance
    virtual syms_t GetIrreps() const;
    virtual void   Display() const;

    size_t GetLMax() const {return itsLMax;}
protected:
    struct NsOnly_t {};
    Atom_EC(int Z, NsOnly_t);       //Populate itsNs only; derived classes build their own occupations.
    struct ValenceOnly_t {};
    //Populate itsNs with the RAW valence shells only (no core, no full-subshell demotion) -- the pseudo-atom
    //(PseudoAtom_EC) config.  \a netCharge != 0 makes an ION: the neutral valence is adjusted by -netCharge
    //electrons (anion adds to the lowest open l, cation removes from the highest occupied l), e.g. F- -> s^2
    //p^6, Na+ -> empty valence.  Unpaired counts are recomputed by Hund's rule on the adjusted valence shell.
    //\a zion is the pseudopotential's q: when it exceeds the neutral valence count the PP keeps SEMICORE shells, which are
    //added as closed shells (Nf) -- the outermost of the noble-gas core, highest n first (Na q9: 2s^2 2p^6; Mn q15: 3s^2 3p^6;
    //Li q3: 1s^2).  A q that is not a whole-shell semicore split THROWS.
    Atom_EC(int Z, int netCharge, int zion, ValenceOnly_t);
    void BuildNROccupations();
    void SetSplitOccupations(sym_t sp, sym_t su, int NCore, int gp, int gu, int Npair, int Nu);
    //! Distribute the atom's \a nup unpaired electrons into itsNs.Nu (Hund's rule on the valence shell,
    //! spilling a single leftover into the next-lower partially-filled l).  Shared by both ctors.
    void AssignUnpaired(int Z, int nup);

    friend class ::ElectronConfigurationTests;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    size_t itsLMax,itsLValance;
    std::map<Irrep,size_t> itsOccupations; //Spin polarized list;
    std::map<Irrep,size_t> itsUnpolOccupations; //Spin un polarized list;
};

//! Pseudo-atom (NR) electron configuration: ONLY the valence shells are occupied, with NO core and NO
//! folding of a filled valence subshell into the core (e.g. Si -> s^2 p^2 over the s and p atomic irreps).
//! This is the configuration a pseudopotential atom wants -- the core electrons are replaced by the
//! pseudopotential, so the SCF carries only the Zion valence electrons over the (core-free) spectrum.
export class PseudoAtom_EC : public Atom_EC
{
public:
    //! \a netCharge is the pseudo-ion's net charge (0 = neutral, -1 = anion e.g. F-, +1 = cation e.g. Na+);
    //! the EC then carries (Zion - netCharge) valence electrons.  \a zion is the pseudopotential's q: a SEMICORE q (Na q9,
    //! Mn q15, Li q3) adds the closed shells it keeps (until 2026-10-31 it was ignored, so a q9 Na atom was solved with ONE
    //! electron: the valence was always the periodic table's neutral one).
    //! \a zion < 0 = the valence-only q (the neutral valence count -- every smallest-q GTH entry).
    PseudoAtom_EC(int Z, int netCharge=0, int zion=-1);
};


} // namespace qchem