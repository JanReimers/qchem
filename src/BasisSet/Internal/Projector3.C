// File: BasisSet/Internal/Projector3.C  The 3-centre projection tensor <ab|c> -- ONE type, several realizations.
//
// V1.1: the merge of the molecular ERI3 (a dense std::vector of <ab|c> matrices) and the reciprocal-space
// G_ERI3 (delta support / matrix-free closures).  The NAME says what both contraction directions do
// (user, 2026-08-16): FORWARD projects the density into the fit space, <rho|c> = Sum_ab D_ab <ab|c> --
// a metric projection, the kernel field says which metric -- and the ADJOINT expands fit-space
// coefficients back to an orbital matrix <i|f|j> = Sum_c f_c <ij|c>.  ("Map3" was rejected: in this
// code region a Map is a keyed container -- ΔG_Map -- and this is an operator.  "ERI3" was rejected as
// inaccurate: the overlap-metric tensor has no "repulsion integral" in it.)
//
// A basis PRODUCES this from its own functions (D-free, cached by BasisSetID); a density CONTRACTS its
// matrix D against it, so D never enters qcBasisSet.  Which realization is set is a property of the
// PRODUCING BASIS, not of the scalar type T -- a deliberate non-correlation (the R2.13/V1.27 lesson):
//   - dense    : the molecular/atomic finite-basis realization -- one full <ab|c> matrix per fit
//                function (RAM-hungry, O(Norb^2 x Nfit)).
//   - columns  : the plane-wave delta support -- for each fit function Delta-m, just the index pairs
//                (i,j): G_i-G_j = Delta-m (weight 1), maximally sparse; `kernel` is the per-column
//                diagonal metric applied on contraction (4pi/|G_c|^2 Coulomb, or EMPTY = overlap).
//   - apply*   : the matrix-free GPW realization -- closures that collocate/integrate on grids
//                without ever materializing a per-column tensor.
module;
#include <cassert>
#include <complex>    // std::operator*(double,complex) -- std header, not a Blaze dep
#include <functional>
#include <map>
#include <type_traits> // std::type_identity_t (ContractAdjoint's non-deduced field parameter)
#include <vector>
#include <utility>
export module qchem.BasisSet.Internal.Projector3;
export import qchem.Types;
import qchem.Blaze;    // rvec_t, hmat_t<T> + complex/double arithmetic

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

//! \brief The 3-centre projection tensor \f$\langle ab|c\rangle\f$ between an orbital basis (indices a,b)
//! and an auxiliary fit basis (index c) -- the currency of every density-fitted Coulomb/Vxc build.
//!
//! FORWARD (the density's side): project \f$D\f$ into the fit space, \f$\langle\rho|c\rangle=\sum_{ab}
//! D_{ab}\langle ab|c\rangle\f$ -- dense realization contracted by the owner (IrrepCD), G-space
//! realizations by \ref Contract.  ADJOINT (the term/fitter side): expand fit coefficients back to an
//! orbital matrix \f$\langle i|f|j\rangle=\sum_c f_c\langle ij|c\rangle\f$ -- dense by the fitter's loop,
//! G-space by \ref ContractAdjoint.  DENSITY-FREE either way: a basis produces this from its own
//! functions; \f$D\f$ never enters the basis.  Exactly ONE realization group is populated (dense, or
//! columns[+kernel], or the apply closures) -- which one is intrinsic to the producing basis lineage.
template <class T> struct Projector3
{
    //! DENSE realization (molecular / atomic): one \f$\langle ab|c\rangle\f$ matrix per fit function c,
    //! in the fit basis's own function order.  RAM-hungry; empty for the G-space lineages.
    std::vector<smat_t<T>> dense;

    struct Column
    {
        ivec3_t                            dm;     //!< \f$\Delta m = G_i-G_j\f$ (the fit function)
        std::vector<std::pair<int,int>>    pairs;  //!< \f$(i,j):G_i-G_j=\Delta m\f$, row-major (\f$i\f$ outer)
    };
    std::vector<Column> columns;   //!< plane-wave delta support; one per fit function (empty for dense/GPW)
    rvec_t              kernel;    //!< per-column \f$4\pi/|G_c|^2\f$ (Coulomb) or EMPTY (overlap, \f$\equiv 1\f$)
    double              volume=0.0;//!< cell volume \f$\Omega\f$ (the \f$1/\Omega\f$ normalisation)

    //! \brief MATRIX-FREE realization (GPW analytic collocation): applies the density-to-\f$\tilde\rho\f$ (or
    //! \f$\to V_H\f$, kernel folded in) map to \f$D\f$ WITHOUT materializing a per-column tensor -- the basis
    //! sets it to a closure that collocates \f$\rho=\sum_{ij}D_{ij}\chi_i\chi_j\f$ analytically on compact
    //! boxes (multi-grid) then FFTs.  When set it takes priority in \ref Contract.  Type-erased so this
    //! leaf names no GPW/grid type.  The Projector3 stays a static-data "spec of required transfers" with an
    //! optional realization.
    std::function<ΔG_Map(const hmat_t<T>& D)> apply;
    //! MATRIX-FREE BACKWARD -- the ADJOINT of \ref apply: assemble the KS matrix \f$\langle i|f|j\rangle=\sum_k
    //! f(G_k)\,W_k(i,j)\f$ from a fitted G-space field \a f (\f$v_{xc}\f$ or \f$V_H\f$).  For a plane-wave basis
    //! this is the Fourier lookup \f$f(m_i-m_j)\f$; for GPW the analytic integrate-back on the fit grid's
    //! REL_CUTOFF ladder.  Always the OVERLAP-metric adjoint (no \ref kernel -- a Coulomb field already carries
    //! \f$4\pi/G^2\f$ from the forward apply).  Because it is produced by \c Overlap3C/Repulsion3C(fitBasis), it
    //! carries the fit grid -- so the KS matrix is 100% consistent with the density that fit basis collocated.
    //! Empty until a basis sets it (forward-only contexts leave it empty).
    std::function<hmat_t<T>(const std::function<T(const ivec3_t&)>& f)> applyAdjoint;
    //! \brief RAW-RASTER forward (doc/GPWPlan 0.5(f2)): \f$\rho_{DM}(r)\f$ on the tensor's integration
    //! raster -- the D-weighted collocated level densities combined SPECTRALLY (zero-pad, no ball
    //! restriction anywhere), the finest level kept raw.  Because \f$D\f$ is PSD, \f$\rho_{DM}=\phi^T D\phi
    //! \ge 0\f$ pointwise to screening precision -- unlike the ball-projected \ref apply whose Gibbs lobes go
    //! negative on sharp products (the XC \f$\rho>0\f$-guard collapse; the C=8 grid-calibration driver).
    //! Empty unless the producing basis has a real-space collocation engine (GPW sets it; plane waves have no
    //! raw representation and leave it empty -- consumers fall back to \ref apply + RhoOnGrid).
    std::function<rvec_t(const hmat_t<T>& D)> applyRaw;
    //! The EXACT adjoint of \ref applyRaw: a real-space field \f$v(r)\f$ on the SAME integration raster ->
    //! the KS matrix \f$\langle i|v|j\rangle\f$ (per-level spectral box-TRUNCATION -- the transpose of the
    //! zero-pad -- then the analytic per-pair gather), so \f$H_{xc}=\partial E_{xc}[\rho_{DM}]/\partial D\f$
    //! holds to machine precision.  Set/empty together with \ref applyRaw.
    std::function<hmat_t<T>(const rvec_t& v)> applyRawAdjoint;
    //! \brief The SAME forward as \ref applyRaw, for a caller holding \f$D\f$ in FACTORED form
    //! \f$D=LL^\dagger\f$ (a thin \f$n\times r\f$ pivoted-Cholesky / natural-orbital factor):
    //! \f$\rho_g=\sum_m\left|[\Phi L]_{gm}\right|^2\f$, i.e. \f$O(n_{pts}nr)\f$ instead of
    //! \f$O(n_{pts}n^2)\f$.
    //!
    //! One integral, two contractions -- WHICH one is the caller's business, because the factor, its
    //! rank memo and the pays-for-itself test are all density-matrix state (ChargeDensity::FactoredRho).
    //! Empty unless the producing basis realises the raw pair over a value table it can left-multiply.
    std::function<rvec_t(const mat_t<T>& L)> applyRawFactored;
};

//! \brief G-space FORWARD contraction: \f$\tilde\rho(\Delta m)=\frac{k_c}\Omega\sum_{(i,j)\in
//! \text{col}(\Delta m)}D_{ij}\f$, with the per-column metric \f$k_c\f$ (\ref Projector3::kernel;
//! \f$k_c\equiv 1\f$ when empty).  The reciprocal-space mirror of the molecular
//! \f$\langle\rho|c\rangle=\sum_{ab}D_{ab}\langle ab|c\rangle\f$ (which the density runs itself over
//! \ref Projector3::dense -- a ΔG_Map is the wrong currency for it).  Dispatches to the matrix-free
//! \ref Projector3::apply when the basis set one.
template <class T> ΔG_Map Contract(const Projector3<T>& g, const hmat_t<T>& D)
{
    // Matrix-free realization (GPW analytic collocation): the closure applies the whole density->rho-tilde
    // (or ->V_H, kernel folded in) map to D directly -- no per-column tensor.  Takes priority when set.
    if (g.apply) return g.apply(D);
    assert(g.dense.empty() && "Contract: a dense Projector3 is contracted by its owner (IrrepCD), not here");
    ΔG_Map rg;
    for (size_t c=0; c<g.columns.size(); ++c)   // plane-wave delta path: k_c/Omega Sum_{(i,j):dm} D_ij
    {
        dcmplx s(0.0);
        for (const auto& [i,j] : g.columns[c].pairs) s += D(i,j);   // Sum_{(i,j): G_i-G_j=dm} D_ij (row-major)
        rg[g.columns[c].dm] = g.kernel.size()==0 ? s/g.volume : g.kernel[c]*s/g.volume;
    }
    return rg;
}

//! The G-space BACKWARD contraction \f$\langle i|f|j\rangle=\sum_k f(G_k)W_k(i,j)\f$ -- the ADJOINT of \ref
//! Contract, building a KS matrix from a fitted G-space field \a f.  Delegates to the matrix-free \ref
//! Projector3::applyAdjoint the basis set (PW Fourier lookup / GPW integrate-back on the FIT grid the tensor
//! was built over) -- so \c ⟨i|v_xc|j⟩ = ContractAdjoint(orb.Overlap3C(vxc_fit), v_xc_tilde), replacing the
//! grid-less \c MakeOverlap(field) bridge (doc/GPWPlan §0e step 2).
// (the field parameter is type_identity_t so a bare lambda still binds -- T deduces from the tensor alone)
template <class T> hmat_t<T> ContractAdjoint(const Projector3<T>& g,
                                             const std::type_identity_t<std::function<T(const ivec3_t&)>>& f)
{
    assert(g.applyAdjoint && "ContractAdjoint: the basis must realise Projector3::applyAdjoint (field->matrix)");
    return g.applyAdjoint(f);
}

//! Frobenius distance between two DENSE realizations (test oracle).
template <class T> double fnorm(const Projector3<T>& a, const Projector3<T>& b);
template <class T> double relative_fnorm(const Projector3<T>& a, const Projector3<T>& b);

} // namespace qchem
