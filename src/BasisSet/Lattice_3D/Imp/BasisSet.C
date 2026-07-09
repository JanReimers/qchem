// File: BasisSet/Lattice_3D/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
#include <cmath>     // lround (fractional k-point -> integer BZ-grid index)
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
// molecular basis is shared (shared_ptr) across the k-blocks.  Rcut folds in lattice images: crystals need
// them for k-dispersion (Rcut=0 makes every k-block identical -- "molecule in a box"); a large-box molecule
// uses the home cell.
GPW_BasisSet::GPW_BasisSet(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                           double densityEcut, double Rcut, double collRcut)
{
    const ivec3_t N=lat.GetLimits();
    for (const auto& kp : lat.MakeKMesh())
    {
        ivec3_t ik(std::lround(kp.k.x*N.x), std::lround(kp.k.y*N.y), std::lround(kp.k.z*N.z));
        // Build the Bloch irrep WITH its BZ weight kp.weight (exactly as PW_BasisSet above) and use the primary
        // sym_t ctor -- the weight carries the Sum_k w_k so the BZ-summed charge/energy are per-cell, not xNk.
        Insert(new GPW_IBS(lat.GetUnitCell(), Symmetry::BlochFactory(N, ik, kp.weight),
                           mol, densityEcut, Rcut, collRcut));                 // mol shared across k-blocks
    }
}

Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, double Rcut, double collRcut)
{
    return new GPW_BasisSet(lat, std::move(mol), densityEcut, Rcut, collRcut);
}

} //namespace
