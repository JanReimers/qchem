// File: Symmetry.C  Abstract interface for symmetries that do not include spin.
module;
#include <cstddef>
#include <Common/pmstream.h>

export module qchem.Symmetry;

export class Symmetry
    : public virtual PMStreamableObject
{
public:
    virtual ~Symmetry() {};
    virtual size_t SequenceIndex() const=0; //Used for op<
    //! Does not include spin degeneracy which is handled separately
    virtual int GetDegeneracy     () const=0;
    virtual int GetPrincipleOffset() const=0; //Add to principle QN.  For atoms this is just l.
    std::string GetLabel          () const;
};
