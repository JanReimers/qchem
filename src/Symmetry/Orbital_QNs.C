// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
module;
#include <iosfwd>

export module qchem.Symmetry.Orbital;
export import qchem.Symmetry.Irrep;
export import qchem.Types;

namespace qchem {

//! \brief An \c Irrep plus a principal quantum number \a n -- the third and outermost level of the
//! orbital-QN hierarchy (\c Symmetry ⊂ \c Irrep ⊂ \c Orbital_QNs).  Fully identifies one orbital block of an
//! atom/molecule/lattice: (principal \a n, spin \a ms, spatial \a sym).  \c SequenceIndex chains \a n onto
//! the \c Irrep key so orbitals sort in shell order; degeneracy is inherited from \c Irrep unchanged.
export struct Orbital_QNs
    : public Irrep
{
    Orbital_QNs() : Irrep(), n(0) {};
    Orbital_QNs(size_t n, Spin ms,const sym_t& sym);
    Orbital_QNs(size_t n, const Irrep&);
    ~Orbital_QNs();

    virtual size_t  SequenceIndex() const;      //!< \c Irrep key chained with the principal QN \a n
    //size_t GetDegeneracy() const; //inherited from Irrep unchanged

    std::ostream& Write(std::ostream&) const;

    size_t  n;                                  //!< principal quantum number of this symmetry block
private:
    static const size_t n_max;                  //!< principal-QN radix used to combine \a n into \c SequenceIndex
};
} // namespace qchem