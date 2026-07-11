// File: Calculation/ValenceBasisGen.C
//
// qchem::GenerateValenceBasis -- generate a valence Gaussian basis for one element under its GTH
// pseudopotential, straight from an atomic pseudo-atom SCF.  The dream (doc/GPWPlan.md sec 1): a grad
// student points this at an element and gets back a well-conditioned, GPW-usable valence .bsd.
//
// WHY the atom calc is the generator/validator (and canned bases are not):
//   * Canned contracted bases (CP2K SZV/DZVP) are optimised for the ANION (e.g. F-), not the neutral atom,
//     and uncontracting them variationally collapses -- so we do NOT copy them.
//   * The MOLECULAR isolated-atom SCF fails for open shells (a degenerate 2p^5 hole rotates freely -> no
//     convergence).  The ATOMIC solver is spherical with the correct l-occupation, so it is the right
//     validator: run the pseudo-atom in a candidate even-tempered window and check the (stable) energy.
//   * The element-wide accuracy pool (Low/Medium/High) spans CORE-scale exponents (~1e4) that are fine for
//     the analytic atomic solver but fatal for the GPW density grid.  So we bring our OWN valence window.
//
// The generated basis is emitted verbatim as the exponents that were validated (single source of truth), in
// Gaussian94 .bsd form -- ready to drop into BasisSetData/ and to transcribe to a CP2K BASIS_SET.
module;
#include <string>
#include <vector>
#include <utility>
export module qchem.ValenceBasisGen;

import qchem.AtomCalculation;   // AtomCalculation, AtomCalcOptions, AtomType (the spherical pseudo-atom solver)

export namespace qchem
{
    //! A valence Gaussian basis to GENERATE for one element under its GTH pseudopotential.  The \c shells are
    //! per-l exponent lists (l, exponents): the atomic pseudo-atom SCF is run in EXACTLY these shells to
    //! validate them (single source of truth -- validate what you emit), then they are emitted verbatim.  The
    //! atomic solver occupies only the l's the CHARGE STATE fills (\c electrons: F- 8, neutral F 7, Na 1), so
    //! higher-l polarization shells ride along un-validated (as intended).  KEEP EXPONENTS DISJOINT ACROSS l:
    //! the molecular Gaussian94 reader merges same-exponent shells (a flagged bug in PG_Cart::IrrepBasisSet),
    //! so shared s/p exponents silently drop functions on load.
    struct ValenceBasisRecipe
    {
        std::string  element;               //!< element symbol, e.g. "F", "Na"
        int          Zion        = 0;       //!< GTH valence charge (0 => the GTH database default)
        int          electrons   = 0;       //!< valence electrons to occupy when validating (0 => Zion, neutral)
        std::string  functional  = "LDA";   //!< GTH parameter set (the GPW/PP stack is LDA today)
        //! The per-l shells: (l, exponents), one uncontracted primitive each.  Validated AND emitted verbatim.
        //! Exponents DISJOINT across shells (see the class note).  Higher-l polarization shells go here too.
        std::vector<std::pair<int,std::vector<double>>> shells;
        std::string  name;                  //!< reserved for future per-recipe naming (file title set at assembly)
    };

    //! The outcome: the atomic pseudo-atom energy in the generated basis (the validation number) and the
    //! Gaussian94 element BLOCK (` El 0 <shells> ****`).  \c converged reports the density gate; the energy is
    //! stable regardless (the degenerate open-shell density need not settle -- see doc/GPWPlan.md).  Following
    //! the BasisSetData convention, a .bsd file is organised by basis-set TYPE (not element): one file holds
    //! many element blocks -- so callers ASSEMBLE the per-element blocks with AssembleBasisFile().
    struct GeneratedBasis { double energy; bool converged; std::string block; };

    //! Validate the recipe (run the pseudo-atom SCF in exactly \c exponents) and emit its element block.
    GeneratedBasis GenerateValenceBasis(const ValenceBasisRecipe&);

    //! Wrap element \a blocks into a complete Gaussian94 .bsd file (the shared comment/`!`/`BASIS=` header,
    //! then every block).  \a name is the file's descriptive title + BASIS= tag.
    std::string AssembleBasisFile(const std::string& name, const std::vector<std::string>& blocks);

    //! A geometric (even-tempered) window: \a N exponents from \a emin to \a emax inclusive.
    std::vector<double> EvenTemperedWindow(int N, double emin, double emax);
}
