// File: Math/GMap.C  G-space coefficient map keyed by a reciprocal-index difference (delta-G).
//
// The shared currency between a plane-wave density (which produces rho-tilde) and a plane-wave basis
// (which consumes it for the Hartree/XC matrices).  Lives in qcMath (the leaf) so qcBasisSet,
// qcChargeDensity, qcLattice_BS and qcHamiltonian can all name the type without a dependency cycle.
// The name emphasises it is a plain CONTAINER keyed by a reciprocal-index difference dG = Gi - Gj, NOT
// a transform -- the "Fourier" provenance lives in the METHODS that produce it (MakeFourierDensity, etc.).
module;
#include <complex>   // std::operator/(complex<double>, double) -- the cross-module operator-visibility stopgap
#include <map>
#include <vector>
#include <utility>
export module qchem.Math.GMap;
import qchem.Types;   // ivec3_t, dcmplx

namespace qchem {

//! Lexicographic comparator so a reciprocal-index triple \f$\Delta m\f$ can key a map.
export struct IVec3Less
{
    bool operator()(const ivec3_t& a, const ivec3_t& b) const
    { if (a.x!=b.x) return a.x<b.x; if (a.y!=b.y) return a.y<b.y; return a.z<b.z; }
};

//! G-space components (density \f$\tilde\rho\f$ or potential \f$\tilde V\f$) keyed by the reciprocal-index
//! difference \f$\Delta m\f$ (\f$\Delta G = B\,\Delta m\f$).
export using ΔG_Map = std::map<ivec3_t, dcmplx, IVec3Less>;

//! \brief The DENSITY-FREE plane-wave three-centre "integrals" \f$\langle G_i G_j|\Delta m\rangle =
//! \tfrac1\Omega\,\delta_{\Delta m,\,G_i-G_j}\f$ -- the reciprocal-space analogue of the molecular \c ERI3
//! tensor \f$\langle ab|c\rangle\f$, exploiting the Kronecker delta so it is stored SPARSELY.
//!
//! One bucket per fit function \f$\Delta m\f$ (an orthonormal \f$\{G\}\f$ "fit basis" function), listing the
//! delta's nonzero index pairs \f$(i,j):G_i-G_j=\Delta m\f$ in row-major (\f$i\f$ outer, \f$j\f$ inner) build
//! order.  A plane-wave BASIS produces this from its own \f$\{G\}\f$ (no density); a plane-wave DENSITY then
//! contracts it against its matrix \f$D\f$ (\ref ContractFourierGather) to get \f$\tilde\rho\f$ -- so the
//! density matrix never crosses into the basis (mirrors molecule: basis owns \c ERI3, density owns \f$D\f$).
export struct FourierGather
{
    std::map<ivec3_t, std::vector<std::pair<int,int>>, IVec3Less> support;  //!< \f$\Delta m\to\{(i,j)\}\f$, row-major
    double volume=0.0;                                                      //!< cell volume \f$\Omega\f$ (the \f$1/\Omega\f$ norm)
};

//! \brief Contract a density matrix against the gather: \f$\tilde\rho(\Delta m)=\tfrac1\Omega\sum_{(i,j)\in
//! \text{bucket}(\Delta m)}D_{ij}\f$ -- the reciprocal-space mirror of the molecular \f$\langle\rho|c\rangle=
//! \sum_{ab}D_{ab}\langle ab|c\rangle\f$.  \a D is ANY type answering \c D(i,j) -> \c dcmplx (templated so
//! this leaf module needs no matrix library).  Bit-identical to the old fused loop: same row-major per-bucket
//! left-fold, same final division by \f$\Omega\f$.
export template <class DMat>
ΔG_Map ContractFourierGather(const FourierGather& g, const DMat& D)
{
    ΔG_Map rg;
    for (const auto& [dm, pairs] : g.support)
    {
        dcmplx s(0.0);
        for (const auto& [i,j] : pairs) s += D(i,j);   // Sum_{(i,j): G_i-G_j=dm} D_ij (row-major order)
        rg[dm] = s / g.volume;                          // 1/Omega normalisation (divide, as the old loop did)
    }
    return rg;
}

} // namespace qchem