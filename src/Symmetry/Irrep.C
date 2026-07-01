// File: Irrep.C  Combine Symmetry with Spin.
module;
#include <iosfwd>
export module qchem.Symmetry.Irrep;
export import qchem.Symmetry;
export import qchem.Symmetry.Spin;
import qchem.Streamable;

namespace qchem {

//! \brief A spatial \c Symmetry plus a \c Spin -- the second level of the orbital-QN hierarchy
//! (\c Symmetry ⊂ \c Irrep ⊂ \c Orbital_QNs, each adding one quantum number).  Value type: it holds the
//! spin channel \a ms and a polymorphic spatial-symmetry handle \a sym.  Ordered by a combined
//! \c SequenceIndex (spatial index interleaved with spin) so irreps sort/deduplicate; degeneracy folds the
//! spatial degeneracy with the spin degeneracy.
export struct Irrep
    : public virtual Streamable
{
    Irrep() : ms(Spin::None), sym(0) {};
    Irrep(Spin _ms,const sym_t& _sym);
    ~Irrep();

    virtual size_t  SequenceIndex() const;      //!< combined (spatial, spin) ordering key

    friend bool operator<(const Irrep& a, const Irrep& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    virtual size_t GetDegeneracy() const;       //!< spatial degeneracy times the spin degeneracy

    std::ostream& Write(std::ostream&) const;

    Spin  ms;                                   //!< spin channel (Up/Down/None)
    sym_t sym;                                  //!< spatial symmetry (polymorphic \c Symmetry handle)

    static const size_t ms_max;                 //!< spin radix used to combine spin into \c SequenceIndex
};



} // namespace qchem