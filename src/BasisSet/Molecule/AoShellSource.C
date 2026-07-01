// File: BasisSet/Molecule/AoShellSource.C  Capability: an orbital basis that can yield its AO shells.
module;
#include <vector>
export module qchem.BasisSet.Molecule.AoShellSource;
export import qchem.Symmetry.OperationRep;   // Symmetry::AoShell

export namespace qchem::BasisSet::Molecule
{

//! \brief Capability interface for an orbital IBS that can produce its AO shells for point-group SALC
//! adaptation.  \c PG::SymmetryAdapt depends ONLY on this abstraction: it does a single abstract-to-abstract
//! cast and calls \c GetAoShells(), rather than \c dynamic_cast-ing to a concrete evaluator data struct
//! (\c PGData / \c SphData).  So a new orbital-basis delivery becomes symmetry-adaptable by implementing one
//! method -- no per-delivery \c if(cast) in the orchestrator, and no silent \f$\{engine\}\times\{angular\}\f$
//! combinatorics trap (a wrong cast can no longer succeed).  Pure interface, no data.
class AoShellSource
{
public:
    virtual ~AoShellSource() {}
    virtual std::vector<Symmetry::AoShell> GetAoShells() const = 0;
};

} //namespace
