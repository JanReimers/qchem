// File: BasisSet/Molecule/PolarizedGaussian/Imp/SymmetryAdapt.C
module;
#include <memory>
#include <cassert>
module qchem.BasisSet.Molecule.PolarizedGaussian.SymmetryAdapt;
import qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;        // ExtractAoShells, ClusterToSymPoints (+ SALC pipeline)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData; // PGData

namespace BasisSet::Molecule::PolarizedGaussian
{
using namespace ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis, const Cluster& cl, double tol)
{
    // The molecular orbital basis is a single PG Orbital_IBS (IS-A PGData).
    const ::BasisSet::Real_OIBS* rawIBS = nullptr;
    const PGData*                pg     = nullptr;
    for (auto ibs : rawBasis->Iterate<const ::BasisSet::Real_OIBS>())
    {
        pg = dynamic_cast<const PGData*>(ibs);
        if (pg) { rawIBS = ibs; break; }
    }
    assert(rawIBS && pg && "PG SymmetryAdapt: no PolarizedGaussian orbital IBS in the basis");

    auto shells = ExtractAoShells(*pg);
    auto pts    = ClusterToSymPoints(cl);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new ::BasisSet::Molecule::SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
