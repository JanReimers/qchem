// File: Pseudopotential/PseudopotentialSet.C  A NAMED pseudopotential choice: for each element, which valence variant (q) of the
// GTH database the run uses, for which functional.
//
// WHY IT EXISTS (doc/CalcNetworkReview.md §4, option D).  The choice of pseudopotential is made ONCE per run, by the run setup, and
// every consumer -- the Hamiltonian's PP terms, the seed library, the basis-provenance check, the XC-grid sharpness gather, the +U
// projector, the electron count, the restart fingerprint -- must agree with it.  A value handed to all of them guarantees that; each
// consumer re-reading gth_potentials.json on its own does not.  This module is the one reader of the choice (data:
// Data/pseudopotentials.json); the potentials themselves stay GetGTH's business.
module;
#include <string>
#include <utility>
#include <vector>
export module qchem.Pseudopotential.Set;

import qchem.Structure;   // Structure -- the atoms a run's vocabulary is read off

export namespace qchem::Pseudopotential
{

class PseudopotentialSet
{
public:
    //! The named set of \c pseudopotentials.json; THROWS listing the known names on a miss.
    static PseudopotentialSet Get(const std::string& name);
    //! Every set name in the data file, in file order.
    static std::vector<std::string> Names();

    const std::string& Name()       const {return itsName;}
    const std::string& Functional() const {return itsFunctional;}   //!< the XC functional the entries were generated for ("LDA", ...)

    //! Does the set resolve \a element: an explicit entry, or (the RULE) any GTH entry for the set's functional.
    bool Covers(const std::string& element) const;
    //! q, the valence electrons kept per atom of \a element (a key of GetGTH): the set's explicit entry if it has one, else its RULE, the
    //! SMALLEST q the GTH database holds for \a element at the set's functional (the valence-only, largest-core variant).  So every element
    //! the database covers resolves; the explicit entries are the exceptions to the rule (OpenWork §3 "A deck cannot name Te with a PP",
    //! 2026-11-11: the old per-element list was the rule written out by hand, and every element it lacked was refused).  THROWS naming the
    //! set when neither covers \a element.
    int  Valence(const std::string& element) const;
    //! The XC functional \a element's resolved potential was FITTED with ("PADE" for the GTH LDA entries; GTH_PP::fitted).  A query for the
    //! policy that pairs or warns on a run whose functional differs (pending a user ruling; nothing reads it yet).  THROWS as \c Valence.
    std::string FittedFunctional(const std::string& element) const;
    //! This set with \a overrides applied: each `(element, q)` replaces the set's q for that element, or ADDS the element when the set lacks it
    //! (a stated pseudopotential valence, e.g. Mn q15 semicore for one run).  The set's name and functional are kept.
    PseudopotentialSet With(const std::vector<std::pair<std::string,int>>& overrides) const;
    //! `(element, q)` for each of \a elements, in that order -- the \c SolidCalcOptions::species form.  THROWS on an uncovered element.
    std::vector<std::pair<std::string,int>> Species(const std::vector<std::string>& elements) const;
    //! `(element, q)` for the distinct elements of \a st's atoms, in order of first appearance -- what the PP terms, the sharpness
    //! gather and the electron count of a run over \a st need (never the whole set).  THROWS on an uncovered element.
    std::vector<std::pair<std::string,int>> Species(const Structure& st) const;

private:
    PseudopotentialSet(std::string name, std::string functional, std::vector<std::pair<std::string,int>> valence)
        : itsName(std::move(name)), itsFunctional(std::move(functional)), itsValence(std::move(valence)) {}
    friend struct SetReader;   // the reader in the implementation unit builds one
    std::string                             itsName, itsFunctional;
    std::vector<std::pair<std::string,int>> itsValence;
};

} // namespace qchem::Pseudopotential
