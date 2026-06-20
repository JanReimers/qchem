// File: BasisSet/Molecule/PolarizedGaussian1/Imp/SymmetryAdapt.C
module;
#include <memory>
#include <cassert>
module qchem.BasisSet.Molecule.PolarizedGaussian1.SymmetryAdapt;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Symmetry;        // ExtractAoShells, ClusterToSymPoints (+ SALC pipeline)
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.PGData; // PGData

namespace BasisSet::Molecule::PolarizedGaussian1
{

::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis, const Cluster& cl, double tol)
{
    // The molecular orbital basis is a single PG1 Orbital_IBS (IS-A PGData).
    const ::BasisSet::Real_OIBS* rawIBS = nullptr;
    const PGData*                pg     = nullptr;
    for (auto ibs : rawBasis->Iterate<const ::BasisSet::Real_OIBS>())
    {
        pg = dynamic_cast<const PGData*>(ibs);
        if (pg) { rawIBS = ibs; break; }
    }
    assert(rawIBS && pg && "PG1 SymmetryAdapt: no PolarizedGaussian1 orbital IBS in the basis");

    auto shells = ExtractAoShells(*pg);
    auto pts    = ClusterToSymPoints(cl);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new ::BasisSet::Molecule::SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
