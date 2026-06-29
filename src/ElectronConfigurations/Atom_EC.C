// File: ElectronConfigurations/Atom_EC.C  Electron configuration for atoms.
module;
#include <map>
#include "forward.H"

export module qchem.ElectronConfiguration.AtomNR;
export import qchem.ElectronConfiguration;

import qchem.ElectronConfiguration.ElectronCounts;
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
    Atom_EC(int Z, ValenceOnly_t);  //Populate itsNs with the RAW valence shells only (no core, no
                                    //full-subshell demotion) -- the pseudo-atom (PseudoAtom_EC) config.
    void BuildNROccupations();
    void SetSplitOccupations(sym_t sp, sym_t su, int NCore, int gp, int gu, int Npair, int Nu);
    //! Distribute the atom's \a nup unpaired electrons into itsNs.Nu (Hund's rule on the valence shell,
    //! spilling a single leftover into the next-lower partially-filled l).  Shared by both ctors.
    void AssignUnpaired(int Z, int nup);

    friend class ElectronConfigurationTests;

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
    PseudoAtom_EC(int Z);
};

