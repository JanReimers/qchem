// File: BasisSet/Molecule/PG_Cart/Imp/SymmetryAdapt.C
// (Doc + the dataflow diagram live on the SymmetryAdapt declaration in the module interface.)
module;
#include <memory>
#include <vector>
#include <stdexcept>
module qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;
import qchem.BasisSet.Molecule.PG_Cart.Symmetry;        // StructureToSymPoints, Centroid (+ the SALC pipeline)
import qchem.BasisSet.Molecule.AoShellSource;           // AoShellSource::GetAoShells -- the one capability we need

namespace qchem::BasisSet::Molecule::PG_Cart
{

::qchem::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::qchem::BasisSet::BasisSet<double>> rawBasis, const Structure& st, double tol)
{
    // Iterate the AoShellSource face directly: the iterator does the (encapsulated) cast, and because an
    // AoShellSource IS-A Real_OIBS we also get the raw IBS for the decorator with no cast-back.  No client
    // dynamic_cast, no PGData/SphData knowledge, no per-delivery branch.  A delivery that can't honour a
    // correct AO-shell layout (e.g. libcint-spherical, S3b) throws from GetAoShells (see SphericalSALCPlan).
    const ::qchem::BasisSet::Real_OIBS*  rawIBS = nullptr;
    std::vector<Symmetry::AoShell>       shells;
    for (auto src : rawBasis->Iterate<const Molecule::AoShellSource>())
        { rawIBS = src; shells = src->GetAoShells(); break; }   // src IS-A Real_OIBS: plain upcast

    if (!rawIBS)
        throw std::runtime_error("PG::SymmetryAdapt: the orbital basis is not symmetry-adaptable "
                                 "(no AoShellSource IBS in the basis).");

    auto pts    = StructureToSymPoints(st);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new ::qchem::BasisSet::Molecule::SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
