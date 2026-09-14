// File: Spin.C  The spin quantum number -- the innermost QN of the orbital-QN hierarchy.
module;
#include <vector>

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

    //! \brief THE IMPOSED SPIN SUBGROUP (doc/CleanupCandidates.md V1.37; doc/BasisSetTaxonomyPlan.md §1.4).
    //!
    //! Spin is a FACTOR of the symmetry group G until it is not.  "Polarized" and "unpolarized" are not
    //! properties of a wave function or a density -- they name WHICH SUBGROUP of the spin factor a run
    //! holds, the same KIND of decision as imposing a point group (a rung on the SSB descent ladder:
    //! impose -> analyse -> release).  Every wave function and density is ONE composite over full
    //! \c Irrep s (spatial ⊗ \c Spin); the subgroup decides only which spin irreps that composite is built
    //! from, and the degeneracy fold (\c Irrep::GetDegeneracy) is what keeps the unpolarized run cheap.
    //!
    //! | imposed | spin group | irreps | \c Spin |
    //! |---|---|---|---|
    //! | \c UnPolarized | full spin rotation SU(2) | one doublet, degeneracy 2 | \c None |
    //! | \c Polarized   | rotations about z, U(1)_z (collinear) | two 1-D irreps \f$m_s=\pm\tfrac12\f$ | \c Up / \c Down |
    //!
    //! (Non-collinear / spin-orbit -- spin INSIDE G, the double group -- is the row that does not exist yet;
    //! a flat composite over double-group irreps needs no new container, which is why the container types
    //! went.)  A factory / policy ARGUMENT, never a type: every Hamiltonian / wave-function factory takes it.
    enum class SpinGroup {UnPolarized, Polarized};

    //! The spin irreps the subgroup resolves -- the channels a composite over full \c Irrep s is built from,
    //! in the order they are built (and hence the block order of every composite density).
    inline std::vector<Spin> SpinIrreps(SpinGroup g)
    {
        return g==SpinGroup::Polarized ? std::vector<Spin>{Spin::Up, Spin::Down} : std::vector<Spin>{Spin::None};
    }
}
