// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
module;
#include <iosfwd>

export module qchem.Symmetry.Orbital;
export import qchem.Symmetry.Irrep;
export import qchem.Types;

namespace qchem {

//! \brief An \c Irrep plus a principal quantum number \a n -- the third and outermost level of the
//! orbital-QN hierarchy (\c Symmetry ⊂ \c Irrep ⊂ \c Orbital_QNs).  Fully identifies one orbital block of an
//! atom/molecule/lattice: (principal \a n, spin \a ms, spatial \a sym).  Orbitals sort LEXICOGRAPHICALLY, the
//! \c Irrep key first and then \a n, so they sort in shell order within a block; degeneracy is inherited from
//! \c Irrep unchanged.
//!
//! ⚠ \a n is UNBOUNDED (lane H, 2026-11-12).  The key used to be the radix pack Irrep*300+(n-1) behind an assert,
//! so a plane-wave block with more than 300 bands (any PW run at a production cutoff) made level (irrep I, n 301)
//! the SAME key as (irrep I+1, n 1): dead assert in Release, a silent collision in \c EnergyLevels' map, so the
//! eigen-row table paired the wrong levels.  There is no radix any more, and \c SequenceIndex is \c Irrep's (the
//! BLOCK key), which is what every caller comparing blocks wants.
export struct Orbital_QNs
    : public Irrep
{
    Orbital_QNs() : Irrep(), n(0) {};
    Orbital_QNs(size_t n, Spin ms,const sym_t& sym);
    Orbital_QNs(size_t n, const Irrep&);
    ~Orbital_QNs();

    //! (block, n) order.  Preferred over \c Irrep's \c operator< for two \c Orbital_QNs (an exact match), so a
    //! map or set keyed on \c Orbital_QNs sorts on both; through an \c Irrep reference it is the block order.
    friend bool operator<(const Orbital_QNs& a, const Orbital_QNs& b)
    {
        const size_t ia=a.Irrep::SequenceIndex(), ib=b.Irrep::SequenceIndex();
        return ia!=ib ? ia<ib : a.n<b.n;
    }
    //size_t GetDegeneracy() const; //inherited from Irrep unchanged

    std::ostream& Write(std::ostream&) const;

    size_t  n;                                  //!< principal quantum number of this symmetry block
};
} // namespace qchem