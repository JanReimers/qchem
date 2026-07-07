// File: BasisSet/Internal/GMap.C  Reciprocal-space (G-space) currencies for the plane-wave / GPW paths.
//
// Lives next to ERI3 (the molecular 3-centre tensor) in qcBasisSet::Internal: these are the reciprocal-space
// analogues, shared by the plane-wave basis (produces them), the charge density (contracts D against them),
// the fitting layer and the Hamiltonian terms.  It sits in qcBasisSet (not the qcMath leaf) so it can use
// Blaze directly (rvec_t, chmat_t, complex arithmetic) without the cross-module operator-visibility stopgaps.
module;
#include <complex>   // std::operator*(double,complex) / operator/(complex,double) -- std header, not a Blaze dep
#include <map>
#include <vector>
#include <utility>
export module qchem.BasisSet.Internal.GMap;
import qchem.Types;   // ivec3_t, dcmplx
import qchem.Blaze;   // rvec_t, chmat_t + complex/double arithmetic (visible here; the qcMath leaf lacked it)

export namespace qchem {

//! Lexicographic comparator so a reciprocal-index triple \f$\Delta m\f$ can key a map.
struct IVec3Less
{
    bool operator()(const ivec3_t& a, const ivec3_t& b) const
    { if (a.x!=b.x) return a.x<b.x; if (a.y!=b.y) return a.y<b.y; return a.z<b.z; }
};

//! G-space components (density \f$\tilde\rho\f$ or potential \f$\tilde V\f$) keyed by the reciprocal-index
//! difference \f$\Delta m\f$ (\f$\Delta G = B\,\Delta m\f$).
using ΔG_Map = std::map<ivec3_t, dcmplx, IVec3Less>;

//! \brief The reciprocal-space three-centre "integrals" \f$\langle G_i G_j|G_c\rangle\f$ -- the G-space
//! analogue of the molecular \c ERI3 tensor \f$\langle ab|c\rangle\f$.
//!
//! One \ref Column per fit function \f$G_c\f$ (a reciprocal-index difference \f$\Delta m\f$ the fit basis
//! carries); for plane waves the integral is the Kronecker delta \f$\delta_{\Delta m,\,G_i-G_j}/\Omega\f$, so
//! a column just lists the orbital index pairs \f$(i,j):G_i-G_j=\Delta m\f$ (weight 1) -- maximally SPARSE.
//! (For GPW the Gaussian product spreads over many \f$G_c\f$, so a column densifies / is realised by grid
//! collocation; the interface is unchanged.)  The \ref kernel is the per-column diagonal metric applied on
//! contraction: \f$4\pi/|G_c|^2\f$ for the Coulomb (Repulsion3C) tensor, or EMPTY for the overlap (Overlap3C)
//! tensor (kernel \f$\equiv 1\f$).  DENSITY-FREE: a basis produces this from its own \f$\{G\}\f$; a density
//! then contracts its matrix \f$D\f$ against it (\ref ContractG_ERI3), so \f$D\f$ never enters the basis.
struct G_ERI3
{
    struct Column
    {
        ivec3_t                            dm;     //!< \f$\Delta m = G_i-G_j\f$ (the fit function)
        std::vector<std::pair<int,int>>    pairs;  //!< \f$(i,j):G_i-G_j=\Delta m\f$, row-major (\f$i\f$ outer)
    };
    std::vector<Column> columns;   //!< one per fit function; the delta support
    rvec_t              kernel;    //!< per-column \f$4\pi/|G_c|^2\f$ (Coulomb) or EMPTY (overlap, \f$\equiv 1\f$)
    double              volume=0.0;//!< cell volume \f$\Omega\f$ (the \f$1/\Omega\f$ normalisation)
};

//! \brief Contract a density matrix against the gather: \f$\tilde\rho(\Delta m)=\frac{k_c}\Omega\sum_{(i,j)\in
//! \text{col}(\Delta m)}D_{ij}\f$, with the per-column metric \f$k_c\f$ (\ref G_ERI3::kernel; \f$k_c\equiv 1\f$
//! when empty).  The reciprocal-space mirror of the molecular \f$\langle\rho|c\rangle=\sum_{ab}D_{ab}\langle
//! ab|c\rangle\f$.  Bit-identical to the old fused loop when the kernel is empty (same row-major per-column
//! fold, same final division by \f$\Omega\f$).
ΔG_Map ContractG_ERI3(const G_ERI3& g, const chmat_t& D)
{
    ΔG_Map rg;
    for (size_t c=0; c<g.columns.size(); ++c)
    {
        dcmplx s(0.0);
        for (const auto& [i,j] : g.columns[c].pairs) s += D(i,j);   // Sum_{(i,j): G_i-G_j=dm} D_ij (row-major)
        rg[g.columns[c].dm] = g.kernel.size()==0 ? s/g.volume : g.kernel[c]*s/g.volume;
    }
    return rg;
}

} // namespace qchem
