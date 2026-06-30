// File: BasisSet/Lattice_3D/Internal/KPlusG.C  Shared reciprocal vectors K=k+G with angular helpers.
//
// Every plane-wave-based lattice basis set (PlaneWave KB nonlocal, APW/LAPW sphere terms) needs the same
// per-plane-wave geometry: the Cartesian reciprocal vector K=k+G, its magnitude |K|, and the angle
// cos(gamma) between two of them that drives the (2l+1)P_l(cos gamma) Legendre sums.  This bundles that
// geometry behind one class so the cos(gamma) guard/clamp and the small-|K| tolerance live in one place.
module;
#include <algorithm>
#include <vector>

export module qchem.BasisSet.Lattice_3D.Internal.KPlusG;
import qchem.ReciprocalLattice;   // UnitCell (ToCartesian)
import qchem.Vector3D;            // norm, operator* (dot)
import qchem.Types;               // rvec3_t, rvec_t, ivec3_t

export namespace qchem::BasisSet::Lattice_3D::Internal
{
//! |k+G| at or below this is treated as the zero vector (its direction is undefined, so cos(gamma)=1).
inline constexpr double kZeroTol = 1e-12;   // TODO: Make External

//! \brief The set of reciprocal vectors \f$K=k+G\f$ (Cartesian) for a plane-wave basis, with their
//! magnitudes and the inter-vector angle \f$\cos\gamma\f$ used by every angular (Legendre) sum.
//! Construct once per basis set / k-point from the cutoff G-set (see Internal::BuildGs).
class KPlusG
{
public:
    KPlusG(const UnitCell& B, const rvec3_t& k, const std::vector<ivec3_t>& G)
        : itsK(G.size()), itsNorm(G.size())
    {
        for (size_t i=0; i<G.size(); i++) { itsK[i]=B.ToCartesian(k+G[i]); itsNorm[i]=norm(itsK[i]); }
    }
    size_t         size()         const {return itsK.size();}
    const rvec3_t& K   (size_t i) const {return itsK[i];}      //!< Cartesian \f$k+G_i\f$.
    double         Norm(size_t i) const {return itsNorm[i];}   //!< \f$|k+G_i|\f$.

    //! \f$\cos\gamma_{ij}\f$ between \f$K_i\f$ and \f$K_j\f$, guarded (a zero vector has no direction, so
    //! the result is 1) and clamped to \f$[-1,1]\f$ against round-off.
    double CosGamma(size_t i, size_t j) const
    {
        double c=(itsNorm[i]>kZeroTol && itsNorm[j]>kZeroTol) ? (itsK[i]*itsK[j])/(itsNorm[i]*itsNorm[j])
                                                              : 1.0;
        return std::max(-1.0, std::min(1.0, c));
    }
private:
    std::vector<rvec3_t> itsK;     //!< Cartesian k+G per plane wave.
    rvec_t               itsNorm;  //!< |k+G| per plane wave.
};
}
