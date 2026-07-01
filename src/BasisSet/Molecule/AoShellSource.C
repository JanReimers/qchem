// File: BasisSet/Molecule/AoShellSource.C  A real orbital basis that can yield its AO shells for SALC.
module;
#include <vector>
export module qchem.BasisSet.Molecule.AoShellSource;
export import qchem.Symmetry.OperationRep;   // Symmetry::AoShell
export import qchem.BasisSet.Orbital_1E_IBS; // Real_OIBS (the orbital-IBS interface this refines)

export namespace qchem::BasisSet::Molecule
{

//! \brief A real orbital IBS that can additionally produce its AO shells for point-group SALC adaptation.
//! It REFINES \c Real_OIBS (is-a orbital 1E basis) rather than being an orthogonal capability, so a client
//! can \c Iterate<AoShellSource> to get exactly the adaptable bases (the iterator does the cast) and still
//! use each as its \c Real_OIBS with no cast-back.  \c PG::SymmetryAdapt depends only on this: it never
//! \c dynamic_cast-s to a concrete evaluator data struct (\c PGData/\c SphData), so a new orbital-basis
//! delivery becomes symmetry-adaptable by implementing one method -- no per-delivery \c if(cast) in the
//! orchestrator, and no silent \f$\{engine\}\times\{angular\}\f$ combinatorics trap.  Deliveries that cannot
//! honour a correct AO-shell layout (e.g. libcint-spherical, whose convention is unmatched) THROW from
//! \c GetAoShells rather than mislead.  Kept off the templated \c Orbital_1E_IBS<T> base on purpose: AO
//! shells are real-space, meaningless for the \c dcmplx (plane-wave) instantiation.
class AoShellSource : public virtual Real_OIBS
{
public:
    virtual std::vector<Symmetry::AoShell> GetAoShells() const = 0;
};

} //namespace
