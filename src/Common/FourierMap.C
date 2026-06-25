// File: Common/FourierMap.C  G-space coefficient map keyed by a reciprocal-index difference.
//
// The shared currency between a plane-wave density (which produces rho-tilde) and a plane-wave basis
// (which consumes it for the Hartree/XC matrices).  Lives in Common so qcBasisSet, qcChargeDensity,
// qcLattice_BS and qcHamiltonian can all name the type without a dependency cycle.
module;
#include <map>
export module qchem.FourierMap;
import qchem.Types;   // ivec3_t, dcmplx

//! Lexicographic comparator so a reciprocal-index triple \f$\Delta m\f$ can key a map.
export struct IVec3Less
{
    bool operator()(const ivec3_t& a, const ivec3_t& b) const
    { if (a.x!=b.x) return a.x<b.x; if (a.y!=b.y) return a.y<b.y; return a.z<b.z; }
};

//! G-space components (density \f$\tilde\rho\f$ or potential \f$\tilde V\f$) keyed by the reciprocal-index
//! difference \f$\Delta m\f$ (\f$\Delta G = B\,\Delta m\f$).
export using FourierMap = std::map<ivec3_t, dcmplx, IVec3Less>;
