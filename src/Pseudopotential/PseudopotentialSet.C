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

    bool Covers(const std::string& element) const;
    //! q, the valence electrons kept per atom of \a element (a key of GetGTH).  THROWS naming the set when it does not cover \a element.
    int  Valence(const std::string& element) const;
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
