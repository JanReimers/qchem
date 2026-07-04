// File: Symmetry/Molecule/Irrep.C  A 1-D molecular point-group irrep (carries its Mulliken label).
//
// The Symmetry object a SymmetryAdapted_IBS carries: GetLabel() (= Write()) returns the
// Mulliken symbol ("A1", "B2", "Ag", ...), which is what decorates the orbital eigenvalue /
// occupation tables.  Abelian groups only, so the degeneracy is 1.  SequenceIndex orders the
// irreps (their column order in the character table).
module;
#include <string>
#include <iosfwd>
export module qchem.Symmetry.Molecule.Irrep;
export import qchem.Symmetry;

namespace qchem::Symmetry::Molecule {

export class Irrep
    : public virtual qchem::Symmetry::Symmetry
{
public:
    Irrep(const std::string& label, size_t index) : itsLabel(label), itsIndex(index) {}

    virtual size_t SequenceIndex    () const {return itsIndex;}
    virtual size_t GetDegeneracy    () const {return 1;}      // abelian: 1-D irreps
    virtual size_t GetPrincipleOffset() const {return 0;}
    virtual std::ostream& Write(std::ostream&) const;         // prints the Mulliken label

private:
    std::string itsLabel;
    size_t      itsIndex;
};

} // namespace qchem::Symmetry::Molecule