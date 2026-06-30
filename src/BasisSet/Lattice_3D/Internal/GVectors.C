// File: BasisSet/Lattice_3D/Internal/GVectors.C  Shared plane-wave G-vector set from an energy cutoff.
module;
#include <cassert>
#include <cmath>
#include <vector>

export module qchem.BasisSet.Lattice_3D.Internal.GVectors;
import qchem.ReciprocalLattice;   // ReciprocalLattice, UnitCell
import qchem.Types;               // ivec3_t, rvec3_t

export namespace qchem::BasisSet::Lattice_3D::Internal
{
//! The reciprocal-lattice index triples \f$m\f$ (\f$G=B\,m\f$) in the plane-wave cutoff set
//! \f$\{G : \tfrac12|k+G|^2 < E_{cut}\}\f$, shared by every plane-wave-based lattice basis set
//! (PlaneWave / APW / LAPW).  \a k is the fractional crystal momentum.
inline std::vector<ivec3_t> BuildGs(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut)
{
    assert(Ecut>0.0);
    const UnitCell& B=recip.GetCell();
    double Gmax=std::sqrt(2*Ecut)+B.GetDistance(k);     // widen the search sphere by |k|
    std::vector<ivec3_t> Gs;
    for (const ivec3_t& m : recip.GetGVectors(Gmax))
    {
        double kG=B.GetDistance(k+m);                   // |k+G|
        if (0.5*kG*kG < Ecut) Gs.push_back(m);
    }
    return Gs;
}
}
