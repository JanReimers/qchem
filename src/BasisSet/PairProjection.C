// File: BasisSet/PairProjection.C  The orbital-PAIR to fit-basis projection <ab|c> -- ONE type, several realizations.
//
// NAMED FOR WHAT IT IS (V-ONE-OF-N, user 2026-11-07; was Projector3): it projects an orbital-pair product onto a fit
// basis, and by realization it is a 3-index (dense <ab|c>), 2-index (PW columns) or 1-index-per-call (the
// closures) object, so "3" named only its first realization (doc/Records/VOneOfNDesign.md §4).
//
// V1.1: the merge of the molecular ERI3 (a dense std::vector of <ab|c> matrices) and the reciprocal-space
// G_ERI3 (delta support / matrix-free closures).  The NAME says what both contraction directions do
// (user, 2026-08-16): FORWARD projects the density into the fit space, <rho|c> = Sum_ab D_ab <ab|c> --
// an OVERLAP projection (V-ONE-OF-N 3c: no kernel in any forward) -- and the ADJOINT expands fit-space
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
//                (i,j): G_i-G_j = Delta-m (weight 1), maximally sparse.  (The per-column Coulomb `kernel` went in
//                V-ONE-OF-N 3c: a periodic Repulsion3C carries only its adjoint, k(G) applied there.)
//   - apply*   : the matrix-free GPW realization -- closures that collocate/integrate on grids
//                without ever materializing a per-column tensor.
module;
#include <stdexcept>   // the realization-mismatch throws (Contract, ContractAdjoint, the integrators)
#include <cassert>
#include <algorithm>   // std::max (MaxAbs)
#include <cmath>       // std::abs (MaxAbs)
#include <complex>    // std::operator*(double,complex) -- std header, not a Blaze dep
#include <functional>
#include <map>
#include <type_traits> // std::type_identity_t (ContractAdjoint's non-deduced field parameter)
#include <vector>
#include <utility>
// ★ IT IS DELIBERATELY NOT `.Internal.` (promoted 2026-09-10) -- and the label was not merely unhelpful,
// it was ILLEGAL.  CLAUDE.md: *"Only re-export modules that don't have `.Internal.` in the module name; the
// goal is to avoid importing internals across library boundaries."*  Two PUBLIC modules were re-exporting
// this one -- `qchem.BasisSet.Orbital_DFT_IBS` and `qchem.BasisSet.GMap` -- which means `PairProjection<T>` and
// `ΔG_Map` have been de facto public vocabulary for as long as `Overlap3C` / `Repulsion3C` have named them
// in their signatures.  The label was a claim the interface had already falsified.
//
// ⇒ Its real consumers are in OTHER families -- `qcFitting` (the delta fitter holds one per block),
// `qcChargeDensity` (`FourierDensity`), `qcHamiltonian` (the pair sampler) -- so under the V1.20 ruling
// (`.Internal.` marks the LIBRARY-FAMILY boundary) this was the same defect as V1.20b's `GMap`, and
// promoting `GMap` alone was only half the fix: the promoted module still re-exported this Internal one.
//
// ⚠ THE DESIGN TENSION IS UNCHANGED AND STILL OPEN, exactly as it is for `GMap`: a 3-centre tensor handle
// is an implementation currency that ideally would not appear in high-level abstract faces, but the
// `Overlap3C` family has to return SOMETHING and that something varies with the basis
// (`doc/CleanupCandidates.md` D2 + V1.20b).  The visibility is the deliverable here, not the cure.
export module qchem.BasisSet.PairProjection;
export import qchem.Types;
export import qchem.Projections;  // ΔG_Map, its field algebra, and the one-of-N record (moved there, V-ONE-OF-N 2/4)
export import qchem.Mesh.Integrator;   // MatrixForward<T>/MatrixAdjoint<T> -- the two halves the raw pair realises
import qchem.Blaze;    // rvec_t, hmat_t<T> + complex/double arithmetic

export namespace qchem {

//! \brief A Bloch-\f$q\f$ G-space field -- \f$\delta\tilde\rho(G+q)\f$ or \f$\delta\tilde V(G+q)\f$ of a transition
//! density -- keyed by the INTEGER G index \f$\Delta m\f$, with the one \f$q\f$ (fractional reciprocal coordinates)
//! every key is shifted by (doc/LinearResponsePlan.md §3d, ruling Q7).
//!
//! ★ A DISTINCT TYPE, NOT A BARE \c ΔG_Map THAT "MEANS" G+q (pin 20: ask what the object MEANS).  A q != 0 map
//! added to a ground-state \f$\tilde\rho\f$, or gathered as if it were lattice-periodic, would compile and be
//! wrong; with its own type that mistake is a build error.  No field algebra is offered: nothing adds two of
//! them yet, and two maps at different q must never be added.
struct ΔGq_Map
{
    rvec3_t q{0,0,0};   //!< the wave vector every key is shifted by (fractional reciprocal coordinates)
    ΔG_Map  c;          //!< the coefficient of \f$e^{i(G+q)\cdot r}\f$, keyed by G's integer index
};

//! \brief The 3-centre projection tensor \f$\langle ab|c\rangle\f$ between an orbital basis (indices a,b)
//! and an auxiliary fit basis (index c) -- the currency of every density-fitted Coulomb/Vxc build.
//!
//! FORWARD (the density's side): project \f$D\f$ into the fit space, \f$\langle\rho|c\rangle=\sum_{ab}
//! D_{ab}\langle ab|c\rangle\f$ -- dense realization contracted by the owner (IrrepCD), G-space
//! realizations by \ref Contract.  ADJOINT (the term/fitter side): expand fit coefficients back to an
//! orbital matrix \f$\langle i|f|j\rangle=\sum_c f_c\langle ij|c\rangle\f$ -- dense by the fitter's loop,
//! G-space by \ref ContractAdjoint.  DENSITY-FREE either way: a basis produces this from its own
//! functions; \f$D\f$ never enters the basis.  Exactly ONE realization group is populated (dense, or
//! columns, or the apply closures; a periodic Repulsion3C is adjoint-only) -- which one is intrinsic to the producing basis lineage.
template <class T> struct PairProjection
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
    double              volume=0.0;//!< cell volume \f$\Omega\f$ (the \f$1/\Omega\f$ normalisation)

    //! \brief MATRIX-FREE realization (GPW analytic collocation): applies the density-to-\f$\tilde\rho\f$
    //! map to \f$D\f$ WITHOUT materializing a per-column tensor -- the basis
    //! sets it to a closure that collocates \f$\rho=\sum_{ij}D_{ij}\chi_i\chi_j\f$ analytically on compact
    //! boxes (multi-grid) then FFTs.  When set it takes priority in \ref Contract.  Type-erased so this
    //! leaf names no GPW/grid type.  The PairProjection stays a static-data "spec of required transfers" with an
    //! optional realization.
    std::function<ΔG_Map(const hmat_t<T>& D)> apply;
    //! MATRIX-FREE BACKWARD -- the ADJOINT of \ref apply: assemble the KS matrix \f$\langle i|f|j\rangle=\sum_k
    //! f(G_k)\,W_k(i,j)\f$ from a fitted G-space field \a f (\f$v_{xc}\f$ or \f$V_H\f$).  For a plane-wave basis
    //! this is the Fourier lookup \f$f(m_i-m_j)\f$; for GPW the analytic integrate-back on the fit grid's
    //! REL_CUTOFF ladder.  On \c Overlap3C the overlap-metric adjoint; on a periodic \c Repulsion3C the COULOMB one,
    //! \f$k(G)\langle i|G|j\rangle\f$ (the Hartree term's physics round, V-ONE-OF-N 3b).  Because it is produced by \c Overlap3C/Repulsion3C(fitBasis), it
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

//! \brief G-space FORWARD contraction: \f$\tilde\rho(\Delta m)=\frac1\Omega\sum_{(i,j)\in
//! \text{col}(\Delta m)}D_{ij}\f$ (an overlap projection: no kernel, V-ONE-OF-N 3c).  The reciprocal-space mirror of the molecular
//! \f$\langle\rho|c\rangle=\sum_{ab}D_{ab}\langle ab|c\rangle\f$ (which the density runs itself over
//! \ref PairProjection::dense -- a ΔG_Map is the wrong currency for it).  Dispatches to the matrix-free
//! \ref PairProjection::apply when the basis set one.
template <class T> ΔG_Map Contract(const PairProjection<T>& g, const hmat_t<T>& D)
{
    // Matrix-free realization (GPW analytic collocation): the closure applies the whole density->rho-tilde
    // map to D directly -- no per-column tensor.  Takes priority when set.
    if (g.apply) return g.apply(D);
    // A MISMATCH THROWS (V-ONE-OF-N): an assert is dead in Release, and a dense projection contracted here would return an
    // empty map, i.e. a silent zero density.  Once per call, so it costs nothing.
    if (!g.dense.empty())
        throw std::logic_error("Contract: this PairProjection is the dense <ab|c> realization; its forward is "
                               "DensePairProjectionIntegrator::Forward, not the G-space Contract");
    // An adjoint-only tensor (a periodic Repulsion3C since V-ONE-OF-N 3c: the kernel is the physics round's) has no forward.
    if (g.columns.empty())
        throw std::logic_error("Contract: this PairProjection realises no G-space forward (no columns, no apply): "
                               "a density projects through Overlap3C, not Repulsion3C");
    ΔG_Map rg;
    for (size_t c=0; c<g.columns.size(); ++c)   // plane-wave delta path: 1/Omega Sum_{(i,j):dm} D_ij
    {
        dcmplx s(0.0);
        for (const auto& [i,j] : g.columns[c].pairs) s += D(i,j);   // Sum_{(i,j): G_i-G_j=dm} D_ij (row-major)
        rg[g.columns[c].dm] = s/g.volume;
    }
    return rg;
}

//! The G-space BACKWARD contraction \f$\langle i|f|j\rangle=\sum_k f(G_k)W_k(i,j)\f$ -- the ADJOINT of \ref
//! Contract, building a KS matrix from a fitted G-space field \a f.  Delegates to the matrix-free \ref
//! PairProjection::applyAdjoint the basis set (PW Fourier lookup / GPW integrate-back on the FIT grid the tensor
//! was built over) -- so \c ⟨i|v_xc|j⟩ = ContractAdjoint(orb.Overlap3C(vxc_fit), v_xc_tilde), replacing the
//! grid-less \c MakeOverlap(field) bridge (doc/GPWPlan §0e step 2).
// (the field parameter is type_identity_t so a bare lambda still binds -- T deduces from the tensor alone)
template <class T> hmat_t<T> ContractAdjoint(const PairProjection<T>& g,
                                             const std::type_identity_t<std::function<T(const ivec3_t&)>>& f)
{
    if (!g.applyAdjoint)
        throw std::logic_error("ContractAdjoint: this PairProjection does not realise the G-space adjoint (applyAdjoint); "
                               "a dense one is DensePairProjectionIntegrator::Adjoint");
    return g.applyAdjoint(f);
}

//! Frobenius distance between two DENSE realizations (test oracle).
template <class T> double fnorm(const PairProjection<T>& a, const PairProjection<T>& b);
template <class T> double relative_fnorm(const PairProjection<T>& a, const PairProjection<T>& b);


//! \brief The SCREENED realization of \c qcMesh::MatrixForward + \c qcMesh::MatrixAdjoint: a
//! \c PairProjection's raw forward/adjoint pair, presented behind the mesh-level faces.
//!
//! ⚠ IT DERIVES FROM THE TWO HALVES, NOT FROM A PAIR FACE (2026-09-09).  There was a
//! \c qcMesh::MatrixIntegrator<T> naming the pair; it was deleted because nothing held it (user: *"nobody
//! should need both sides … it just adds confusion"*).  The no-mismatch guarantee is that ONE of these
//! objects is constructed and its two halves are handed to the two clients that each need one -- a
//! property of construction, not of a type.  Hence the ctor check below, which is the only place both
//! directions are looked at together.
//!
//! ★ THE INTERFACE FIT WAS NOT ENGINEERED, IT WAS DISCOVERED (2026-09-08).  The two faces were designed at
//! the mesh level from the adjointness requirement alone, and their signatures turn out to be EXACTLY the
//! ones \c applyRaw and \c applyRawAdjoint already had:
//! \code
//!     std::function<rvec_t   (const hmat_t<T>& D)> applyRaw;         // Forward(D) -> rho
//!     std::function<hmat_t<T>(const rvec_t&    v)> applyRawAdjoint;  // Adjoint(v) -> matrix
//! \endcode
//! So this is a THIN ADAPTER, not a reimplementation -- which is the strongest evidence available that the
//! face is the right one: two independently-derived descriptions of the same seam agreed.
//!
//! WHAT IT TRUNCATES, and why that is the whole point.  Unlike the dense sibling this route screens per
//! pair and box-truncates per multigrid level -- but it does so IDENTICALLY in both directions, which is
//! what makes \f$H_{xc}=\partial E_{xc}[\rho_{DM}]/\partial D\f$ exact for the TRUNCATED functional.  The
//! dense realization is its natural reference: same answer to quadrature accuracy, at
//! \f$O(n_{pts}n^2)\f$ instead of the screened pair count.
//!
//! \warning It borrows the \c PairProjection -- which is owned by the basis's integral cache and outlives any
//! assembly -- so this object is a VIEW and must not outlive it.
template <class T> class ScreenedMatrixIntegrator
    : public virtual qcMesh::MatrixForward<T>
    , public virtual qcMesh::MatrixAdjoint<T>
{
public:
    //! \a g must realise the raw pair; \a weights are the integration raster's (they size the coefficient
    //! axis, and supply the default \c Integrate rule \f$\sum_a w_a f_a\f$).
    //!
    //! ★ \a rule OVERRIDES THAT QUADRATURE, AND IT EXISTS FOR A PINNED SUMMATION ORDER (2026-09-09).  The
    //! GPW raster's own energy rule is \f$(\sum_a f_a)\,\Omega/N\f$ -- algebraically the same as uniform
    //! weights, DIFFERENT in the last bits, and the periodic energies are pinned to it at ten digits.  A
    //! client that must integrate the raster's way therefore supplies the rule instead of accepting ours.
    //! ⚠ WHY A CONSTRUCTOR ARGUMENT AND NOT A SECOND `Integrate`: the integration rule is a DECISION about
    //! which discrete functional the run minimises, so it belongs to whoever constructs -- the one actor
    //! that knows.  Leaving it to the call site is how an object ends up with two defensible answers and a
    //! caller picking by accident, which is the whole failure this family exists to prevent.
    //!
    //! THROWS when the raw pair is absent rather than assembling half a route: a \c PairProjection with a
    //! forward and no adjoint cannot be an integrator, and finding that out at the first \c Adjoint call --
    //! inside an SCF iteration -- is strictly worse than finding it out here.
    ScreenedMatrixIntegrator(const PairProjection<T>& g, rvec_t weights,
                             std::function<double(const rvec_t&)> rule={})
        : itsG(g), itsW(std::move(weights)), itsRule(std::move(rule))
    {
        if (!g.applyRaw || !g.applyRawAdjoint)
            throw std::runtime_error("ScreenedMatrixIntegrator: this PairProjection does not realise the raw "
                "forward/adjoint pair.  An integrator is the PAIR -- a tensor that can collocate but not "
                "integrate back (or the reverse) is not one, and pretending otherwise is exactly the "
                "mismatch this interface exists to prevent.");
    }

    //! Un-hide the base's FACTORED overload before overriding this one (see qchem.Mesh.Integrator's
    //! warning): without it `Forward(L)` converts a thin factor into an hmat_t adaptor and Blaze throws.
    using qcMesh::MatrixForward<T>::Forward;
    virtual Projections Forward(const hmat_t<T>& D) const override {return Projections(itsG.applyRaw(D));}
    //! The FAST factored path when the tensor has one; otherwise the base's default (form \f$LL^\dagger\f$,
    //! delegate) answers -- correct either way, which is what makes it not a capability.
    virtual Projections Forward(const mat_t<T>& L) const override
    {
        if (itsG.applyRawFactored) return Projections(itsG.applyRawFactored(L));
        return qcMesh::MatrixForward<T>::Forward(L);
    }
    virtual hmat_t<T> Adjoint(const rvec_t& v)    const override {return itsG.applyRawAdjoint(v);}
    virtual double    Integrate(const rvec_t& f)  const override
    {
        assert(f.size()==itsW.size() && "ScreenedMatrixIntegrator::Integrate: one value per raster point");
        if (itsRule) return itsRule(f);                  // the client's pinned rule (see the ctor)
        double s=0.0;
        for (size_t g=0; g<f.size(); ++g) s+=itsW[g]*f[g];
        return s;
    }
    virtual size_t NumCoefficients() const override {return itsW.size();}

private:
    const PairProjection<T>& itsG;   //!< borrowed: owned by the basis's integral cache
    rvec_t               itsW;   //!< the integration raster's weights (also the coefficient-axis length)
    //! The constructing client's own energy quadrature, when it has a pinned one; empty => use \c itsW.
    std::function<double(const rvec_t&)> itsRule;
};

//! \brief The DENSE, ANALYTIC realization of the forward/adjoint pair -- a \c PairProjection whose
//! \ref PairProjection::dense holds one full \f$\langle ab|c\rangle\f$ matrix per fit function (the
//! molecular Gaussian auxiliary basis; Coulomb or overlap metric, the tensor says which).
//!
//! ★ THIS IS R1.0q's FIRST ANALYTIC INSTANCE (2026-09-19), AND R1.0p's PROOF BY CONSTRUCTION: nothing here
//! is a grid.  The coefficient axis is the fit basis's FUNCTIONS, the forward is
//! \f$\rho_a=\sum_{ab}D_{ab}\langle ab|a\rangle\f$ and the adjoint \f$\sum_a v_a\langle ij|a\rangle\f$ --
//! the same tensor contracted the other way, so \f$H=\partial E/\partial D\f$ holds exactly, with no
//! screening to keep consistent because there is none.  The two loops are the ones \c FiniteIrrepCD::
//! GetRepulsion3C and \c ConstrainedFF::Repulsion used to run on the SAME borrowed tensor from opposite
//! sides of a library boundary, agreeing by convention; they now come off ONE object, and the summation
//! order is kept so the move is bit-identical.
//!
//! \warning A VIEW: it borrows \a g, which the producing basis's integral cache owns for the run.
template <class T> class DensePairProjectionIntegrator
    : public virtual qcMesh::MatrixForward<T>
    , public virtual qcMesh::MatrixAdjoint<T>
{
public:
    //! \a integrals are the fit functions' own \f$\langle f_a|1\rangle\f$ (\c Charge()), which is what
    //! \c Integrate contracts a coefficient vector against -- function integrals, never point weights.
    DensePairProjectionIntegrator(const PairProjection<T>& g, rvec_t integrals)
        : itsG(g), itsInt(std::move(integrals))
    {
        if (itsG.dense.empty())
            throw std::runtime_error("DensePairProjectionIntegrator: this PairProjection has no dense realization -- "
                "it is the analytic (molecular/atomic) integrator and needs one <ab|c> matrix per fit function.");
        if (itsInt.size()!=itsG.dense.size())
            throw std::runtime_error("DensePairProjectionIntegrator: one function integral per fit function");
    }

    using qcMesh::MatrixForward<T>::Forward;   // un-hide the factored overload (see qchem.Mesh.Integrator)
    //! \f$\rho_a=\sum_{ab}D_{ab}\langle ab|a\rangle\f$ -- the loop \c FiniteIrrepCD::GetRepulsion3C ran.
    virtual Projections Forward(const hmat_t<T>& D) const override
    {
        rvec_t ret(itsG.dense.size());
        for (size_t a=0; a<itsG.dense.size(); ++a)
            ret[a]=IReal(blazem::sum(D % itsG.dense[a]));
        return Projections(std::move(ret));
    }
    //! \f$\langle i|v|j\rangle=\sum_a v_a\langle ij|a\rangle\f$ -- the loop \c ConstrainedFF::Repulsion ran.
    virtual hmat_t<T> Adjoint(const rvec_t& v) const override
    {
        assert(v.size()==itsG.dense.size() && "DensePairProjectionIntegrator::Adjoint: one coefficient per fit function");
        hmat_t<T> J=blazem::zeroH<T>(itsG.dense.empty() ? 0 : itsG.dense[0].rows());
        for (size_t a=0; a<v.size(); ++a) J+=v[a]*itsG.dense[a];
        return J;
    }
    virtual double Integrate(const rvec_t& f) const override
    {
        assert(f.size()==itsInt.size() && "DensePairProjectionIntegrator::Integrate: one value per fit function");
        double s=0.0;
        for (size_t a=0; a<f.size(); ++a) s+=itsInt[a]*f[a];
        return s;
    }
    virtual size_t NumCoefficients() const override {return itsG.dense.size();}

private:
    template <class U> static double IReal(const U& x)
    {
        if constexpr (std::is_floating_point_v<U>) return x; else return std::real(x);
    }
    const PairProjection<T>& itsG;     //!< borrowed: owned by the basis's integral cache
    rvec_t               itsInt;   //!< \f$\langle f_a|1\rangle\f$ per fit function
};

//! \brief The G-SPACE realization: a \c PairProjection whose forward is \ref Contract (plane-wave columns, or GPW's
//! matrix-free \ref PairProjection::apply), and whose adjoint is \ref ContractAdjoint (\ref PairProjection::applyAdjoint).
//! The forward returns the \c ΔG_Map realization of \c Projections: \f$\tilde\rho(G)\f$ through an \c Overlap3C.
//!
//! ★ IT IMPLEMENTS \c MatrixForward ONLY.  Its adjoint takes the plane-wave fitter's OWN coefficients, a \c ΔG_Map, so it
//! is a typed member here and not \c MatrixAdjoint (whose input is a vector): fit coefficients are each fitter's private
//! type (user ruling 2026-11-08 (b), doc/Records/VOneOfN_D2Design.md §2a), and nothing generic calls a G-space adjoint.
//!
//! \warning A VIEW: it borrows \a g, which the producing basis's integral cache owns for the run.
//! \a U is the orbital BLOCK's scalar; the fit axis is complex either way (a real TRIM block contracts against the run's
//! complex fit axis, 3c-3), so the tensor is always a \c PairProjection<dcmplx>.
template <class U> class GSpacePairProjectionIntegrator
    : public virtual qcMesh::MatrixForward<U>
{
public:
    //! THROWS unless \a g realises both G-space directions (a forward -- columns or \c apply -- and \c
    //! applyAdjoint): an integrator is the PAIR, and finding a missing half at the first SCF call is worse.
    explicit GSpacePairProjectionIntegrator(const PairProjection<dcmplx>& g) : itsG(g)
    {
        if (!g.dense.empty() || !g.applyAdjoint || (!g.apply && g.columns.empty()))
            throw std::runtime_error("GSpacePairProjectionIntegrator: this PairProjection does not realise the G-space "
                "forward (columns or apply) and adjoint (applyAdjoint) pair");
    }

    using qcMesh::MatrixForward<U>::Forward;   // un-hide the factored overload (see qchem.Mesh.Integrator)
    virtual Projections Forward(const hmat_t<U>& D) const override {return Projections(Contract(itsG, D));}
    //! No native factored path: a rebuilt \f$LL^\dagger\f$ would only cost a matrix and roundoff (see MatrixForward::FactorPays).
    virtual bool FactorPays() const override {return false;}
    //! \f$\langle i|f|j\rangle=\sum_G f(G)\,W_G(i,j)\f$ for a field \a f on the same \f$\{G\}\f$; a G absent from the map reads
    //! as 0.  Complex for either block scalar (the caller narrows a real TRIM block, as every periodic term does).
    hmat_t<dcmplx> Adjoint(const ΔG_Map& f) const
    {
        return ContractAdjoint(itsG, [&f](const ivec3_t& dm)->dcmplx
            { auto it=f.find(dm); return it==f.end() ? dcmplx(0.0) : it->second; });
    }
    //! The number of plane-wave columns.  THROWS for the matrix-free realization: its G set is the closure's
    //! own business and has no length to report (a ΔG_Map is sized by its keys, not by a count).
    virtual size_t NumCoefficients() const override
    {
        if (itsG.columns.empty())
            throw std::logic_error("GSpacePairProjectionIntegrator::NumCoefficients: a matrix-free G-space projection has no fixed coefficient count");
        return itsG.columns.size();
    }

private:
    const PairProjection<dcmplx>& itsG;   //!< borrowed: owned by the basis's integral cache
};

} // namespace qchem