// File: BasisSet/Molecule/PolarizedGaussian/Imp/SymmetryAdaptedBasisSet.C
module;
#include <string>
#include <memory>
#include <cassert>
module qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;
import qchem.Symmetry.MolecularIrrep;     // MolecularIrrep (the labelled irrep symmetry)
import qchem.Blaze;                        // submatrix
import qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;          // ExtractAoShells, ClusterToSymPoints
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;   // PGData

namespace BasisSet::Molecule
{
using namespace ::BasisSet::Molecule::PolarizedGaussian;            // ExtractAoShells, ClusterToSymPoints

SymmetryAdaptedBasisSet::SymmetryAdaptedBasisSet(const ::BasisSet::Orbital_1E_IBS<double>* raw,
                                                 const Symmetry::SALCs& salc)
{
    size_t nAO = salc.O.rows();
    for (size_t r=0; r+1<salc.blockStart.size(); ++r)        // one non-empty irrep block -> one IBS
    {
        size_t start = salc.blockStart[r], dG = salc.blockStart[r+1]-start;
        if (dG==0) continue;                                 // skip irreps absent from this basis
        rmat_t            Or    = blazem::submatrix(salc.O, 0, start, nAO, dG);
        const std::string label = salc.irrep[start];         // all columns of the block share it
        sym_t             sym(new MolecularIrrep(label, r));
        this->Insert(new ::BasisSet::SymmetryAdapted_IBS(raw, Or, label, sym));
    }
}

SymmetryAdaptedBasisSet* SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis,
                                       const Cluster& cl, double tol)
{
    // The molecular orbital basis is a single PG Orbital_IBS (IS-A PGData).
    const ::BasisSet::Real_OIBS* rawIBS = nullptr;
    const PGData*                pg     = nullptr;
    for (auto ibs : rawBasis->Iterate<const ::BasisSet::Real_OIBS>())
    {
        pg = dynamic_cast<const PGData*>(ibs);
        if (pg) { rawIBS = ibs; break; }
    }
    assert(rawIBS && pg && "SymmetryAdapt: no PolarizedGaussian orbital IBS in the basis");

    auto shells = ExtractAoShells(*pg);
    auto pts    = ClusterToSymPoints(cl);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
