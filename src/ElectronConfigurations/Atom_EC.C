// File: ElectronConfigurations/Atom_EC.C  Electron configuration for atoms.
module;
#include <map>
#include <string>
#include <vector>
#include "forward.H"

export module qchem.ElectronConfiguration.AtomNR;
export import qchem.ElectronConfiguration;
export import qchem.Outcome;   // ParseConfiguration: a configuration TEXT is user input

import qchem.ElectronConfiguration.ElectronCounts;

namespace qchem {
const int Nshell=8;

//! One subshell of a configuration: principal number \c n, angular momentum \c l, and its \c count of electrons.
export struct ShellOccupation { int n, l, count; };
export using ShellOccupations = std::vector<ShellOccupation>;

//! \brief An electron configuration as STATED in NIST's notation ("[Ar] 3d^6 4s^1"; the "^" is optional): the noble-gas core \c [X]
//! expanded, and the subshells written after it (2026-11-14, OpenWork row "All-electron ION seeds": an ion is solved in a STATED
//! configuration, since Fe+ is 3d^6 4s^1 where its isoelectronic neutral Mn is 3d^5 4s^2).
export struct StatedConfiguration
{
    std::string      coreSymbol;   //!< the [X] noble gas ("Ar"), empty when none is written
    ShellOccupations core;         //!< the [X] subshells, expanded
    ShellOccupations named;        //!< the subshells written after it, as written
    int CoreElectrons () const;
    int NamedElectrons() const;
    ShellOccupations All() const;  //!< core then named
};
//! Parse \a text (NIST's notation).  Refused (an Outcome: the text is user input): an unknown core, a token that is not n l count, l past f,
//! n <= l, a count past the subshell's capacity 2(2l+1), or a subshell written twice.
export Outcome<StatedConfiguration, std::string> ParseConfiguration(const std::string& text);
//! \brief Why the radial atom cannot hold \a shells, or empty when it can.  The radial solver fills each l block from its LOWEST level up,
//! so in each l only the highest occupied n may be open, every lower occupied n is full, and the occupied n's are consecutive (no hole
//! below an occupied level); \a fromLowestN (all-electron) also requires each l to start at n = l+1 (a pseudo-atom's l blocks start at its
//! valence).  At most one open subshell per l follows.
export std::string WhyNotRepresentable(const ShellOccupations& shells, bool fromLowestN);
//! NIST's notation for \a shells, sorted by (n, l) ("1s^2 2s^2 ... 3d^6 4s^1").
export std::string ConfigurationText(const ShellOccupations& shells);

export class Atom_EC : public virtual ElectronConfiguration
{
public:
    Atom_EC(int Z);
    //! \brief The atom in a STATED configuration \a shells -- every subshell the SCF carries (all-electron: the whole atom; a pseudo-atom:
    //! the shells beyond its pseudopotential's core).  Closed subshells are core, each l's open one is valence, with Hund's maximal spin
    //! (min(k, 2(2l+1) - k) unpaired; NIST's LSD atoms and ions are high-spin).  \a shells must be representable (\c WhyNotRepresentable
    //! with fromLowestN false): a caller that did not check is a broken invariant, and THROWS std::logic_error.
    explicit Atom_EC(const ShellOccupations& shells);
    //! The configuration this EC occupies, as subshells: each l's \c N[l] electrons filled from n = l+1 + (that l's shells in the core
    //! the EC does not carry: none for an all-electron atom, the pseudopotential's for a pseudo-atom).  What the run report prints.
    ShellOccupations Shells() const;
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
    //Li q3: 1s^2).  When it is SMALLER (rule B, 2026-11-15) the valence's closed subshells go into the core, f^14 first then d^10
    //(I q7: 5s^2 5p^5 over [Kr] 4d^10).  A q reached by neither whole-shell route THROWS.
    Atom_EC(int Z, int netCharge, int zion, ValenceOnly_t);
    void BuildNROccupations();
    void SetSplitOccupations(sym_t sp, sym_t su, int NCore, int gp, int gu, int Npair, int Nu);
    //! Distribute the atom's \a nup unpaired electrons into itsNs.Nu (Hund's rule on the valence shell,
    //! spilling a single leftover into the next-lower partially-filled l).  Shared by both ctors.
    void AssignUnpaired(int Z, int nup);

    friend class ::ElectronConfigurationTests;

    static const int FullShells[Nshell][LMax+2];
    ElCounts itsNs; //Total,core, valance and unpaired counts.
    //! Per l, how many whole subshells lie BELOW the ones itsNs counts: 0 for an all-electron atom; for a pseudo-atom, that l's shells in
    //! the core its pseudopotential replaces (Kr q8: s 3, p 2, d 1 -- the core [Ar] 3d^10 is no prefix of the Madelung filling).
    int itsCoreShells[LMax+1]={0,0,0,0};
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