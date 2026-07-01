// File: BasisSet/Molecule/PG_Cart/Imp/SymmetryAdapt.C
module;
#include <memory>
#include <vector>
#include <stdexcept>
module qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;
import qchem.BasisSet.Molecule.PG_Cart.Symmetry;        // ExtractAoShells(PGData), StructureToSymPoints (+ SALC)
import qchem.BasisSet.Molecule.PG_Spherical.Symmetry;   // ExtractAoShells(SphData) -- the in-house spherical path
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;          // PGData
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;            // SphData

namespace qchem::BasisSet::Molecule::PG_Cart
{
using namespace ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD
namespace Sph = ::qchem::BasisSet::Molecule::Evaluators::PG_Spherical_MnD;

::qchem::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::qchem::BasisSet::BasisSet<double>> rawBasis, const Structure& st, double tol)
{
    // Pick the AoShell extractor by orbital-IBS angular kind: the Cartesian PolarizedGaussian (IS-A PGData)
    // or the in-house spherical basis (IS-A SphData).  NOTE (CLAUDE.md dynamic_cast survey): these cast an
    // abstract Real_OIBS to a CONCRETE data struct -- the established pattern for the data-carrying IBS
    // types (the Cartesian PGData cast predates this); flag for the system-wide cast pass.
    const ::qchem::BasisSet::Real_OIBS*  rawIBS = nullptr;
    std::vector<Symmetry::AoShell>       shells;
    for (auto ibs : rawBasis->Iterate<const ::qchem::BasisSet::Real_OIBS>())
    {
        if (auto* pg = dynamic_cast<const PGData*>(ibs))
            { rawIBS = ibs; shells = ExtractAoShells(*pg); break; }
        if (auto* sph = dynamic_cast<const Sph::SphData*>(ibs))
            { rawIBS = ibs; shells = PG_Spherical::ExtractAoShells(*sph); break; }
    }
    // Guard (not an assert): the libcint-spherical delivery is still unsupported (S3b) -- and it presents
    // AS a PGData with spherical components, a trap the facade guards against upstream (engine=MnD only for
    // symmetry+spherical).  See doc/SphericalSALCPlan.md.
    if (!rawIBS)
        throw std::runtime_error("PG::SymmetryAdapt: no supported orbital IBS -- symmetry adaptation covers "
                                 "the Cartesian PolarizedGaussian (PGData) and in-house spherical (SphData) "
                                 "bases; libcint-spherical is not yet wired.");

    auto pts    = StructureToSymPoints(st);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new ::qchem::BasisSet::Molecule::SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
