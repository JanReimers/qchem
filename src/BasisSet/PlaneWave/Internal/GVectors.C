// File: BasisSet/PlaneWave/Internal/GVectors.C  Shared plane-wave G-vector set from an energy cutoff.
module;
#include <cassert>
#include <cmath>
#include <vector>

export module qchem.BasisSet.PlaneWave.Internal.GVectors;
import qchem.ReciprocalLattice;   // ReciprocalLattice, UnitCell
import qchem.LatticePoints;       // LatticePointsInSphere (on the sizing box)
import qchem.Types;               // ivec3_t, rvec3_t

export namespace qchem::BasisSet::PlaneWave::Internal
{
//! The reciprocal-lattice index triples \f$m\f$ (\f$G=B\,m\f$) in the plane-wave cutoff set
//! \f$\{G : \tfrac12|k+G|^2 < E_{cut}\}\f$, shared by every plane-wave-based lattice basis set
//! (PlaneWave / APW / LAPW).  \a k is the fractional crystal momentum.
//! { m : ½|k + B m|² < Ecut }, k fractional.  The ball is chosen on the reciprocal box's SIZING metric (\c CellGeometry::SizingGeometry;
//! the box itself unless a frozen rule set one), so a strained cell keeps the reference's index set.
inline std::vector<ivec3_t> BuildGs(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut)
{
    assert(Ecut>0.0);
    const CellGeometry B=recip.GetCell().SizingGeometry();
    double Gmax=std::sqrt(2*Ecut)+B.GetDistance(k);     // widen the search sphere by |k|
    std::vector<ivec3_t> Gs;
    for (const ivec3_t& m : LatticePointsInSphere(B, Gmax))
    {
        double kG=B.GetDistance(k+m);                   // |k+G|
        if (0.5*kG*kG < Ecut) Gs.push_back(m);
    }
    return Gs;
}
}
