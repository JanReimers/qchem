// File: BasisSet/Molecule/PG_Cart/Imp/SymmetryAdapt.C
module;
#include <memory>
#include <stdexcept>
module qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;
import qchem.BasisSet.Molecule.PG_Cart.Symmetry;        // ExtractAoShells, StructureToSymPoints (+ SALC pipeline)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData; // PGData

namespace BasisSet::Molecule::PG_Cart
{
using namespace ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis, const Structure& st, double tol)
{
    // The molecular orbital basis is a single PG Orbital_IBS (IS-A PGData).
    const ::BasisSet::Real_OIBS* rawIBS = nullptr;
    const PGData*                pg     = nullptr;
    for (auto ibs : rawBasis->Iterate<const ::BasisSet::Real_OIBS>())
    {
        pg = dynamic_cast<const PGData*>(ibs);
        if (pg) { rawIBS = ibs; break; }
    }
    // Guard (not an assert): symmetry adaptation needs the Cartesian PolarizedGaussian (PGData) shells.
    // The spherical / libcint-spherical deliveries are not yet supported -- see doc/SphericalSALCPlan.md.
    if (!(rawIBS && pg))
        throw std::runtime_error("PG::SymmetryAdapt: basis has no Cartesian PolarizedGaussian (PGData) "
                                 "orbital IBS; symmetry adaptation is only supported for the default "
                                 "Cartesian PG basis (not spherical/libcint).");

    auto shells = ExtractAoShells(*pg);
    auto pts    = StructureToSymPoints(st);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new ::BasisSet::Molecule::SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
