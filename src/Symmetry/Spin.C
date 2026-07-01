// File: Spin.C  The spin quantum number -- the innermost QN of the orbital-QN hierarchy.
module;

export module qchem.Symmetry.Spin;
export import qchem.Types;

export namespace qchem
{
    //! \brief The spin channel of an orbital/irrep.  \c None is the spin-unpolarized case (both channels
    //! collapsed, degeneracy 2); \c Up / \c Down are the polarized channels (degeneracy 1 each).  This is
    //! the innermost quantum number: a spatial \c Symmetry gains spin to become an \c Irrep.
    enum class Spin {Down,None,Up};
    inline bool   IsPolarized  (Spin s) {return !(s==Spin::None);}          //!< true unless \c None
    inline int    GetDegeneracy(Spin s) {return IsPolarized(s) ? 1 : 2;}    //!< 1 per polarized channel, 2 for \c None
    inline size_t SequenceIndex(Spin s) {return static_cast<size_t>(s);}    //!< ordering key (Down<None<Up)
}