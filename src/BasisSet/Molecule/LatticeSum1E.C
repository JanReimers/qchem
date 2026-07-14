// File: BasisSet/Molecule/LatticeSum1E.C  Periodic (lattice-summed) 1-electron integrals of a molecular
// Gaussian basis -- the GPW seam.
//
// A molecular Gaussian orbital basis is a set of contracted Gaussians standing at the atoms.  Placed in a
// periodic cell, the Bloch orbital is the lattice sum  chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R), and its
// one-electron matrices are the corresponding lattice sums of the ordinary (finite) two-centre integrals:
//   S_ij = Sum_R e^{ik.R} <chi_i | chi_j(.-R)> ,  and likewise <p^2> and the nuclear attraction.
// The R=0 term IS the finite integral, so with R={0} these reproduce the molecule's own Overlap/Kinetic/
// Nuclear matrices exactly.
//
// This is the ONE high-level operation GPW needs from the molecular basis: it hands the basis a list of
// real-space lattice translations (Cartesian, MUST include {0,0,0}) and gets back the summed matrices.  The
// Gaussian machinery -- the primitives, the M&D kernels, GaussianRF::AtCenter that places a radial at each
// image -- stays ENCAPSULATED on the molecular side (the caller never sees a radial or an exponent); GPW
// stays purely lattice-side.  A concrete molecular Gaussian basis realises this face and GPW reaches it by an
// abstract->abstract cross-cast of the orbital block.
//
// ENGINE-NEUTRAL SEAM (the integral-engine switch point).  This face is deliberately engine-agnostic: it
// names no integral engine, so it can be realised by ANY molecular Gaussian basis.  Today the McMurchie-
// Davidson Cartesian basis (PG_Cart::Orbital_IBS, Engine::MnD) implements it -- its AtCenter + analytic 2C
// kernels make the per-image sum exact and trivial, ideal for the correctness increment.  The FASTER libCint
// engine (PG_LibCint) can realise the same face (per image, the cross-centre integrals between the home and
// R-shifted shells) -- a perf follow-up.  Because GPW consumes ONLY this abstract face (it takes any
// Molecule::Real_BS and cross-casts, throwing cleanly if unsupported), switching engines is then just the
// Engine argument to Molecule::Factory where the basis is built -- NO change to GPW itself.
//
// GENERAL k (this increment): the phase e^{ik.R} enters as a per-image complex weight, so the sums are
// Hermitian (chmat_t) -- at Gamma every phase is 1 and the result is real (imaginary part exactly zero).  The
// caller supplies the (Rs, phases) pair (phases[r] is the coefficient for Rs[r]; the origin's phase is 1); the
// k-convention (2 pi k_frac . n) stays entirely on the lattice/GPW side, so this face is geometry-only.  A
// physically rigorous periodic nuclear attraction (Ewald, vs. this large-cell single-image limit) is a later
// increment.
//
// FUTURE: (Rs, phases) is a weighted point set -- {R} points + {e^{ik.R}} complex weights -- so once qcMesh
// grows a complex-weight Mesh (cMesh = Mesh<dcmplx>) the pair collapses to a single `const cMesh&` argument
// (Rs -> m.Points(), phases -> m.Weights()).  Kept as an adjacent pair for now to make that swap mechanical.
module;
#include <vector>
export module qchem.BasisSet.Molecule.LatticeSum1E;
import qchem.Structure;   // Structure (the nuclear-attraction centres)
import qchem.Types;       // rvec3_t, cvec_t, chmat_t
import qchem.Blaze;       // matrix machinery

export namespace qchem::BasisSet::Molecule
{

//! \brief The periodic (lattice-summed) 1-electron capability of a molecular Gaussian basis: the general-k
//! Bloch matrices \f$M_{ij}(k)=\sum_R e^{ik\cdot R}\langle\chi_i|\,\hat O\,|\chi_j(\cdot-R)\rangle\f$ over a
//! caller-supplied set of Cartesian lattice translations \a Rs (which MUST contain the origin) and the
//! matching per-image phases \a phases (\c phases[r] weights \c Rs[r]; the origin's phase is 1).  With
//! \c Rs={0} (phase 1) each returns the finite molecule's own matrix; growing \a Rs folds in the periodic
//! images.  Hermitian (\c chmat_t); real at \f$\Gamma\f$ (all phases 1 \f$\Rightarrow\f$ zero imaginary part).
//!
//! \c Make prefix / not self-caching: like the other 1E builders, these are the single-shot producers behind
//! the cached \c Integrals_{Overlap,Kinetic,Nuclear} accessors (\c Overlap()/\c Kinetic()/\c Nuclear() on the
//! GPW orbital IBS cache them in \c DB_Cache keyed on \c BasisSetID, so a static-across-SCF matrix is built
//! once and reused).  The \a Rs/\a phases arguments already mark these as the lattice (periodic) forms.
//!
//! \note \a Rs and \a phases are an adjacent (points, weights) pair -- see the FUTURE note in the file header:
//! they collapse to one \c const \c cMesh& once qcMesh grows a complex-weight \c Mesh.
class LatticeSum1E
{
public:
    virtual ~LatticeSum1E() = default;

    //! \f$S_{ij}=\sum_R e^{ik\cdot R}\langle\chi_i|\chi_j(\cdot-R)\rangle\f$ (normalised).
    virtual chmat_t MakeOverlap(const std::vector<rvec3_t>& Rs, const cvec_t& phases) const = 0;
    //! \f$\langle p^2\rangle_{ij}=\sum_R e^{ik\cdot R}\langle\chi_i|-\nabla^2|\chi_j(\cdot-R)\rangle\f$ (NO 1/2 --
    //! the Hamiltonian applies it, matching \c Integrals_Kinetic).
    virtual chmat_t MakeKinetic(const std::vector<rvec3_t>& Rs, const cvec_t& phases) const = 0;
    //! \f$\sum_R e^{ik\cdot R}\langle\chi_i|\sum_c -Z_c/|r-R_c|\,|\chi_j(\cdot-R)\rangle\f$ over the atoms of \a cl.
    virtual chmat_t MakeNuclear(const std::vector<rvec3_t>& Rs, const cvec_t& phases, const Structure* cl) const = 0;

    //! The largest primitive Gaussian exponent \f$\alpha_{\max}\f$ in the basis -- a scalar RESOLUTION summary
    //! (NOT primitive exposure; the radials stay encapsulated).  A GPW density grid must resolve the sharpest
    //! density feature, the product of the two tightest primitives (a Gaussian of exponent \f$2\alpha_{\max}\f$),
    //! so its minimum plane-wave cutoff is \f$\propto\alpha_{\max}\f$.  This lets the grid consumer floor its own
    //! \c densityEcut from the basis instead of leaving that (basis-dependent) burden on the caller.
    virtual double MaxExponent() const = 0;

    //! \brief The GPW COLLOCATION-ADJOINT (integrate-back), molecular-side and PATCHED -- the KS block
    //! \f$M_{ij}=w\sum_{p\in\mathrm{supp}(i,j)}\overline{\chi_i^k(r_p)}\,V(r_p)\,\chi_j^k(r_p)\f$ over the GPW
    //! density grid \a gridPts, restricted per orbital to the grid points where the Bloch orbital
    //! \f$\chi_i^k=\sum_R e^{ik\cdot R}\chi_i(\cdot-R)\f$ is above the (internal) magnitude-screening tolerance
    //! -- so a pair contracts only on the intersection of its two orbitals' supports (the tight/localized part).
    //! \a Rs / \a phases are the Bloch image set (as in \c MakeOverlap); \a V the real-space potential sampled
    //! on \a gridPts; \a w the uniform quadrature weight \f$\Omega/N_{pts}\f$.  Hermitian (real diagonal), real
    //! at \f$\Gamma\f$.  The Gaussian primitives (centres, exponents) stay ENCAPSULATED -- the screening / support
    //! determination is internal; only the grid + potential cross the seam.  Bit-consistent with the dense grid
    //! contraction \f$w\,\Phi^H(V\!\odot\!\Phi)\f$ to the screening tolerance.  This is the SINGLE-GRID scaffold
    //! for the multi-grid rewrite (doc/GPWPlan.md \S0): per-orbital patches on one grid are ~break-even with the
    //! dense GEMM, but they carry directly to per-exponent grid levels where each pair's patch is small.
    virtual chmat_t MakePotentialMatrix(const rvec3vec_t& gridPts,
                                        const std::vector<rvec3_t>& Rs, const cvec_t& phases,
                                        const rvec_t& V, double w) const = 0;

    //! The coarsest primitive Gaussian exponent \f$\alpha_{\min}\f$ -- the diffuse end (mirrors \c MaxExponent).
    //! Sets the coarsest useful GPW density-grid LEVEL for the multi-grid integrate-back (\c MakePotentialMatrixMG):
    //! the ladder runs from the fine grid (\f$\propto\alpha_{\max}\f$) down to \f$\propto\alpha_{\min}\f$.
    virtual double MinExponent() const = 0;

    //! \brief MULTI-GRID collocation adjoint (integrate-back): each orbital PAIR is contracted on the COARSEST
    //! grid LEVEL that resolves its product exponent \f$\alpha_i+\alpha_j\f$, so diffuse pairs live on small
    //! coarse grids instead of the fine grid dictated by the tightest primitive
    //! (\f$O(n^2N_{pts}^{fine})\to O(\sum_L\mathrm{pairs}(L)\,N_{pts}(L))\f$ -- the GPW gap-closer, CP2K's
    //! \c REL_CUTOFF multigrid).  Levels arrive FINEST-FIRST: \a gridPts_L[L] the level's grid points,
    //! \a V_L[L] the potential restricted to that level, \a w_L[L]\f$=\Omega/N_{pts}(L)\f$ the quadrature weight,
    //! \a ecut_L[L] the level cutoff (descending -- drives the internal pair\f$\to\f$level assignment from the
    //! encapsulated exponents).  \a Rs / \a phases the Bloch image set.  APPROXIMATE vs the single fine grid --
    //! a coarse-level pair carries that level's (adequate-by-construction) grid error; converges to the fine
    //! result as the ladder refines.  With \f$K=1\f$ it reduces to \c MakePotentialMatrix on the fine grid.
    virtual chmat_t MakePotentialMatrixMG(const std::vector<rvec3vec_t>& gridPts_L, const std::vector<double>& ecut_L,
                                          const std::vector<rvec3_t>& Rs, const cvec_t& phases,
                                          const std::vector<rvec_t>& V_L, const std::vector<double>& w_L) const = 0;
};

} //namespace
