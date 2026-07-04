// File: Symmetry/Unit.C  Abstract interface for no symmetry.
module;
#include <iosfwd>
export module qchem.Symmetry.Unit;
export import qchem.Symmetry;
export import qchem.Types;

namespace qchem {

export class UnitQN
    : public virtual Symmetry::Symmetry
{
public:
    UnitQN();
    virtual size_t SequenceIndex() const {return 1;} //Used for op<
    virtual size_t GetDegeneracy() const;
    virtual size_t GetPrincipleOffset() const {return 0;}

    virtual std::ostream& Write(std::ostream&) const;
};



} // namespace qchem