// File: BasisSet/Lattice_3D/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
#include <cmath>     // lround (fractional k-point -> integer BZ-grid index)
module qchem.BasisSet.Lattice_3D.BasisSet;
import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> (the generic list-of-IBS container)
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

} //namespace
