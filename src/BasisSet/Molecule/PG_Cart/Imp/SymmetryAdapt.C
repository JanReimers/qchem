// File: BasisSet/Molecule/PG_Cart/Imp/SymmetryAdapt.C
//! \file
//! \brief Symmetry-adapt a raw molecular basis: \c Structure + \c RawBasisSet \f$\rightarrow\f$ \c SALC_IBS.
//! \image html salc_call_flow.svg "Structure + RawBasisSet -> SALC_IBS dataflow" width=640
//! (source: doc/diagrams/salc_call_flow.svg; add doc/diagrams to Doxyfile IMAGE_PATH to render.)
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
    // The orbital IBS produces its own AO shells (in its own angular convention) through the AoShellSource
    // capability -- ONE abstract-to-abstract cast, no knowledge of the concrete PGData/SphData data structs
    // and no per-delivery branch.  A basis that can't be symmetry-adapted simply isn't an AoShellSource, so
    // the guard below reports it cleanly (e.g. libcint-spherical, S3b; see doc/SphericalSALCPlan.md).
    const ::qchem::BasisSet::Real_OIBS*  rawIBS = nullptr;
    std::vector<Symmetry::AoShell>       shells;
    for (auto ibs : rawBasis->Iterate<const ::qchem::BasisSet::Real_OIBS>())
        if (auto* src = dynamic_cast<const Molecule::AoShellSource*>(ibs))
            { rawIBS = ibs; shells = src->GetAoShells(); break; }

    if (!rawIBS)
        throw std::runtime_error("PG::SymmetryAdapt: the orbital basis is not symmetry-adaptable (no "
                                 "AoShellSource) -- covered today are the Cartesian PolarizedGaussian and "
                                 "in-house spherical bases; libcint-spherical is not yet wired.");

    auto pts    = StructureToSymPoints(st);
    auto grp    = Symmetry::BuildAbelianGroup(pts, tol);
    auto salc   = Symmetry::BuildSALCs(shells, grp, Symmetry::Centroid(pts), tol);

    auto* sab = new ::qchem::BasisSet::Molecule::SymmetryAdaptedBasisSet(rawIBS, salc);
    sab->KeepAlive(rawBasis);                    // self-contained: holds the raw basis alive
    return sab;
}

} //namespace
