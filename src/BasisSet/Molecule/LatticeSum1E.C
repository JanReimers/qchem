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
#include <functional>   // cellphase_t (the Bloch phase of an integer cell offset -- k stays lattice-side)
#include <vector>
export module qchem.BasisSet.Molecule.LatticeSum1E;
import qchem.Structure;      // Structure (the nuclear-attraction centres)
import qchem.UnitCell;       // UnitCell (the collocation grid<->cell map: CollocateDensity)
import qchem.Types;          // rvec3_t, cvec_t, chmat_t, rmat_t, ivec3_t
import qchem.Blaze;          // matrix machinery
import qchem.Math.Angular;   // Math::CartTerm (the Cartesian-monomial expansion of a GaussianFunction)

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

    //! \brief A primitive Cartesian-Gaussian scalar function,
    //! \f$g(r)=\sum_t c_t\,(r-C)^{m_t}\,e^{-\alpha|r-C|^2}\f$ (one shared exponent, a finite Cartesian-
    //! monomial polynomial \f$(r-C)^{m}=x^{m_x}y^{m_y}z^{m_z}\f$ about the centre \f$C\f$) -- the family of
    //! scalar functions a Gaussian basis can integrate ANALYTICALLY.  Pure function-language: the caller may
    //! think of \f$g\f$ as a projector, a moment, anything -- this face only sees a Gaussian.
    struct GaussianFunction
    {
        rvec3_t center;                         //!< \f$C\f$
        double  alpha;                          //!< the shared exponent \f$\alpha\f$
        std::vector<Math::CartTerm> terms;      //!< \f$\{(m_t,c_t)\}\f$ -- the polynomial about \f$C\f$
    };

    //! \brief The lattice-summed overlap of every basis function with ONE Gaussian function \a g:
    //! \f$b_i=\sum_R \mathrm{phases}[R]\,\langle\chi_i|\,g(\cdot-R)\,\rangle\f$ -- the vector (\f$n\f$)
    //! analogue of the matrix \c MakeOverlap, with \a g standing in the \f$\chi_j\f$ slot.  Same
    //! \c (Rs,phases) weighted point set, same magnitude screen (\f$g\f$'s reach from \f$\alpha\f$).
    //! \f$\chi_i\f$ carries its usual normalisation; \a g is integrated RAW (its scale is in \c terms).
    virtual cvec_t MakeOverlap(const std::vector<rvec3_t>& Rs, const cvec_t& phases,
                               const GaussianFunction& g) const = 0;
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

    //! The coarsest primitive Gaussian exponent \f$\alpha_{\min}\f$ -- the diffuse end (mirrors \c MaxExponent).
    //! Sets the coarsest useful GPW density-grid LEVEL for the multi-grid collocation ladder: the levels run
    //! from the fine grid (\f$\propto\alpha_{\max}\f$) down to \f$\propto\alpha_{\min}\f$.
    virtual double MinExponent() const = 0;

    //! \brief The Bloch phase of an INTEGER cell offset \f$n\f$: \f$e^{ik\cdot R_n}\f$.  The k-CONVENTION stays
    //! entirely on the lattice/GPW side (the caller supplies this closure); the molecular side only ever asks
    //! "what is the phase of offset \f$n\f$" for the cross-cell pair offsets it enumerates internally (it cannot
    //! receive a pre-built \c (Rs,phases) pair because the offsets are magnitude-screened per pair, not a fixed
    //! caller-visible set).  \f$\Gamma\f$ = the constant 1.
    using cellphase_t = std::function<dcmplx(const ivec3_t& n)>;

    //! \brief ANALYTIC per-pair density collocation on a MULTI-GRID ladder (CP2K GPW / Quickstep): the grid
    //! density \f$\rho(r)=\sum_{ij}\sum_{R}\,\mathrm{Re}[D_{ij}e^{-ik\cdot R}]\,\chi_i(r)\chi_j(r-R)\f$ --
    //! the cell-periodic density of the Bloch orbital products \f$\chi_i^k\overline{\chi_j^k}\f$ -- with each
    //! (pair, cross-cell offset \f$R\f$) term evaluated ANALYTICALLY on its compact exp-tail box and
    //! MODULO-WRAPPED onto the grid.  The offsets are magnitude-screened (no hard \c Rcut, no Gibbs ringing);
    //! the wrap IS the image sum.  Density-matrix driven, so \f$\int\rho=\mathrm{Tr}(DS^k)\f$ (\f$S^k\f$ the
    //! screened-complete Bloch overlap).  MULTI-GRID (CP2K \c REL_CUTOFF): levels arrive FINEST-FIRST --
    //! \a N_L[L] the level's grid divisions of cell \a A (raster \f$r=A(idx/N)\f$), \a ecut_L[L] its cutoff
    //! (descending); each pair is collocated on the COARSEST level that still resolves its product exponent
    //! \f$\alpha_i+\alpha_j\f$ (the internal assignment -- primitives stay encapsulated), so diffuse pairs live
    //! on small coarse grids: \f$O(n^2N_{pts}^{fine})\to O(\sum_L\mathrm{pairs}(L)N_{pts}(L))\f$.  Returns one
    //! grid density per level (the caller FFTs each and combines \f$\tilde\rho\f$ nested in G-space).  Because
    //! the collocation is analytic (never a sampled orbital), a diffuse pair on its matched coarse grid is
    //! accurate -- the sampling multigrid's fatal aliasing is absent by construction.  \f$K=1\f$ = single grid.
    virtual std::vector<rvec_t> CollocateDensity(const chmat_t& D, const cellphase_t& phase, const UnitCell& A,
                                                 const std::vector<ivec3_t>& N_L,
                                                 const std::vector<double>& ecut_L) const = 0;

    //! \brief The collocation ADJOINT (integrate-back): the KS block \f$h_{ij}=\langle\chi_i^k|V|\chi_j^k\rangle
    //! =\sum_R e^{ik\cdot R}\,w_L\sum_{\text{box}}\chi_i(r)\chi_j(r-R)V(r)\f$, per (pair, offset) on the SAME
    //! compact exp-tail box + modulo-wrap + level assignment as \c CollocateDensity -- so it is the EXACT
    //! adjoint (variational: \f$\int\rho V=\mathrm{Tr}(Dh)\f$ to machine precision).  \a V_L[L] is the real-space
    //! potential on level \f$L\f$'s grid (the caller restricts \f$\tilde V\f$ to each level's \f$\{G\}\f$ --
    //! a SPECTRAL low-pass, no ringing); \f$w_L=\Omega/N_{pts}(L)\f$ is internal.  Only \f$V\f$ is sampled
    //! (weighted by the analytic Gaussians), never the sharp orbital product -> accurate on the matched coarse
    //! grid.  Hermitian; real at \f$\Gamma\f$.
    //! \a relCutoffScale STIFFENS the pair->level requirement for SHARP fields (default 1 = the smooth-field
    //! calibration).  The energy error of a (pair x field) product decays with the SUM of both spectra: for the
    //! smooth \f$V_H/V_{xc}\f$ the field decays like the density, but the LOCAL PSEUDOPOTENTIAL is spectrally
    //! BROAD, so its integrate-back must place each pair on a finer level (~6x) -- while the ultra-diffuse
    //! pairs still fall to deep coarse levels (their own spectra kill the field's tail).
    virtual chmat_t IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                                       const std::vector<ivec3_t>& N_L,
                                       const std::vector<double>& ecut_L, double relCutoffScale=1.0) const = 0;
};

} //namespace
