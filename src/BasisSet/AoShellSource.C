// File: BasisSet/AoShellSource.C  "I am built from atom-centred shells": the capability face a consumer
// selecting functions by (site, l) needs -- a Hubbard MANIFOLD (+U, 2026-09-19), a SALC induction.
//
// It was a pure virtual on the molecular Gaussian block alone (Gaussian::Point::Orbital_1E_IBS::GetAoShells,
// for the point-group SALC builder).  A Bloch block built by Bloch-summing that molecular block is built
// from the SAME shells in the SAME order, so it can honestly answer too -- and the +U term, which sees only
// the abstract orbital block, needs the answer on a face it can cross-cast to (abstract to abstract), not
// on a concrete lattice class.  Hence a face of its own, at the level where both lineages can derive it.
//
// ⚠ It says nothing about angular convention: a shell's rep (Cartesian monomials or real solid harmonics)
// carries that.  A consumer that needs a PARTICULAR convention checks the rep, never the count.
module;
#include <vector>
export module qchem.BasisSet.AoShellSource;
export import qchem.Symmetry.Molecule.OperationRep;   // Symmetry::Molecule::AoShell (+ ShellRep::L)

export namespace qchem::BasisSet
{
class AoShellSource
{
public:
    virtual ~AoShellSource() = default;
    //! \brief This block's shells: centre, angular rep (which knows its \f$l\f$), per-component normalisation,
    //! and the OFFSET of each shell's first component in this block's function ordering.  Deliveries that
    //! cannot honour a correct layout THROW (libcint-spherical, S3b).
    virtual std::vector<Symmetry::Molecule::AoShell> GetAoShells() const = 0;
};
}
