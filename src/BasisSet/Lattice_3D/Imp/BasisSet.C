// File: BasisSet/Lattice_3D/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
#include <cmath>     // lround (fractional k-point -> integer BZ-grid index)
#include <iostream>  // std::cout (the run-start GPW grid diagnostic)
#include <memory>    // std::shared_ptr / std::move (the GPW basis owns the molecular Gaussian basis)
module qchem.BasisSet.Lattice_3D.BasisSet;
import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> (the generic list-of-IBS container)
import qchem.BasisSet.Lattice_3D.GPW_IBS;     // GPW_IBS (the periodic-Gaussian block GPW_BasisSet owns)
import qchem.Symmetry.Factory;                // BlochFactory (the Bloch irrep per k)
import qchem.Types;

namespace qchem::BasisSet::Lattice_3D
{

// ONE plane-wave block per Brillouin-zone k-point: the basis ctor is the single place that enumerates k, so
// the framework's per-irrep loop (each k IS a Bloch irrep) becomes the BZ sum Sum_k w_k.  The KMesh carries
// the points + weights (uniform 1/Nk for an unreduced grid; symmetry-reduced points/weights will plug in
// here later).  N=(1,1,1) -> a single Gamma block.
PW_BasisSet::PW_BasisSet(const ::qchem::Lattice_3D& lat, double Ecut)
{
    ReciprocalLattice recip=lat.Reciprocal();
    const ivec3_t N=lat.GetLimits();
    for (const auto& kp : lat.MakeKMesh())
    {
        ivec3_t ik(std::lround(kp.k.x*N.x), std::lround(kp.k.y*N.y), std::lround(kp.k.z*N.z));
        auto* pw=new PlaneWave_IBS(recip, Symmetry::BlochFactory(N, ik, kp.weight), Ecut);
        Insert(pw);                                 // BasisSetImp takes ownership (no PP: the Vpseudo
                                                    // Hamiltonian term owns the pseudopotential model)
    }
}

Complex_BS* Factory(Type type, const ::qchem::Lattice_3D& lat, double Ecut)
{
    assert(type==Type::PW);
    return new PW_BasisSet(lat, Ecut);
}

// GPW: one Bloch block of periodic Gaussians per Brillouin-zone k-point, built from the molecular basis over
// the cell's atoms -- the exact mirror of PW_BasisSet above (the ctor is the sole place that enumerates k, so
// the framework's per-irrep loop becomes the BZ sum Sum_k w_k).  N=(1,1,1) -> a single Gamma block.  The
// molecular basis is shared (shared_ptr) across the k-blocks.  Lattice images are eps-converged series
// summed inside the molecular seam (no radius exists); homeCellOnly is the finite "molecule in a box" MODE
// (image-free by definition -- every k-block identical), used by the finite==lattice gates.
GPW_BasisSet::GPW_BasisSet(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                           double densityEcut, rvec3_t kShift, CellImages images, double cutoffFactor)
{
    const ivec3_t N=lat.GetLimits();
    const GPW_IBS* first=nullptr;
    for (const auto& kp : lat.MakeKMesh(kShift))
    {
        // Recover the INTEGER grid index (undo the shift first, then round) -- lround(kp.k*N) alone is broken
        // for shift=½ (rounds i+½ to the wrong integer).  ik + kShift reconstruct the exact k in BlochFactory.
        ivec3_t ik(std::lround(kp.k.x*N.x-kShift.x), std::lround(kp.k.y*N.y-kShift.y), std::lround(kp.k.z*N.z-kShift.z));
        // Build the Bloch irrep WITH its BZ weight kp.weight (exactly as PW_BasisSet above) and use the primary
        // sym_t ctor -- the weight carries the Sum_k w_k so the BZ-summed charge/energy are per-cell, not xNk.
        auto* b=new GPW_IBS(lat.GetUnitCell(), Symmetry::BlochFactory(N, ik, kp.weight, kShift),
                            mol, densityEcut, images, cutoffFactor);   // mol shared across k-blocks
        if (!first) first=b;
        Insert(b);
    }
    // GRID DIAGNOSTIC (doc/GPWPlan §0e): basis exponents + every stored grid, once per run.  The grids are
    // k-independent (the density grid/ladder are built at k=0 in every block), so the first block speaks for all.
    if (first) first->ReportGrids(std::cout);
}

Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, rvec3_t kShift, CellImages images, double cutoffFactor)
{
    return new GPW_BasisSet(lat, std::move(mol), densityEcut, kShift, images, cutoffFactor);
}

} //namespace
