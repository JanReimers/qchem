// File: BasisSet/Molecule/LatticeSum1E.C  Periodic (lattice-summed) 1-electron integrals of a molecular
// Gaussian basis -- the GPW seam.
//
// A molecular Gaussian orbital basis is a set of contracted Gaussians standing at the atoms.  Placed in a
// periodic cell, the Bloch orbital is the lattice sum  chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R), and its
// one-electron matrices are the corresponding lattice sums of the ordinary (finite) two-centre integrals:
//   S_ij = Sum_R e^{ik.R} <chi_i | chi_j(.-R)> ,  and likewise <p^2> and the nuclear attraction.
//
// THERE IS NO CUT -- in the R direction (user pin, doc/GPWPlan.md).  A lattice sum is a CONVERGENT SERIES,
// summed to eps by the magnitude screen; it is never truncated at a radius.  The ENUMERATION therefore
// lives HERE, per shell pair (the integrand's owner enumerates: which offsets matter is a function of the
// Gaussian tails -- data only this side owns), exactly as CollocateDensity/IntegratePotential already do.
// No radius, no translation list, no weighted point set crosses the interface: the caller hands the cell
// geometry + a Bloch-phase ORACLE (cellphase_t -- "what is the weight of integer offset n") and gets the
// summed matrices back.  The Gaussian machinery -- the primitives, the M&D kernels, GaussianRF::AtCenter
// that places a radial at each image -- stays ENCAPSULATED on the molecular side (the caller never sees a
// radial or an exponent); GPW stays purely lattice-side (it owns k and the phase convention e^{2 pi i k.n}).
// A concrete molecular Gaussian basis realises this face and GPW reaches it by an abstract->abstract
// cross-cast of the orbital block.
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
// GENERAL k: the phase e^{ik.R} enters as a per-offset complex weight through the oracle, so the sums are
// Hermitian (chmat_t) -- at Gamma every phase is 1 and the result is real (imaginary part exactly zero).
// A physically rigorous periodic nuclear attraction (Ewald, vs. this large-cell single-image limit) is a
// later increment.
//
// (The old "(Rs, phases) -> one cMesh" future note is MOOT for this seam -- no weighted point set crosses
// the interface any more, the stronger form of that cleanup.  KMesh + the quadrature meshes still want the
// Mesh<W> templating; see doc/GPWPlan.md section 5.)
module;
#include <functional>   // cellphase_t (the Bloch phase of an integer cell offset -- k stays lattice-side)
#include <vector>
export module qchem.BasisSet.Molecule.LatticeSum1E;
import qchem.Structure;      // Structure (the nuclear-attraction centres)
import qchem.UnitCell;       // UnitCell (the cell geometry the internal enumeration walks)
import qchem.Types;          // rvec3_t, cvec_t, chmat_t, rmat_t, ivec3_t
import qchem.Blaze;          // matrix machinery
import qchem.Math.Angular;   // Math::CartTerm (the Cartesian-monomial expansion of a GaussianFunction)

export namespace qchem::BasisSet::Molecule
{

//! \brief The periodic (lattice-summed) 1-electron capability of a molecular Gaussian basis: the general-k
//! Bloch matrices \f$M_{ij}(k)=\sum_R e^{ik\cdot R}\langle\chi_i|\,\hat O\,|\chi_j(\cdot-R)\rangle\f$.
//! The series is summed to \f$\varepsilon\f$ INTERNALLY, per shell pair (magnitude screening -- THERE IS NO
//! CUT); the caller supplies only the cell \a A and the Bloch-phase oracle \a phase.  Hermitian
//! (\c chmat_t); real at \f$\Gamma\f$ (all phases 1 \f$\Rightarrow\f$ zero imaginary part).
//!
//! \c Make prefix / not self-caching: like the other 1E builders, these are the single-shot producers behind
//! the cached \c Integrals_{Overlap,Kinetic,Nuclear} accessors (\c Overlap()/\c Kinetic()/\c Nuclear() on the
//! GPW orbital IBS cache them in \c DB_Cache keyed on \c BasisSetID, so a static-across-SCF matrix is built
//! once and reused).  The \a (phase,A) arguments mark these as the lattice (periodic) forms.
class LatticeSum1E
{
public:
    virtual ~LatticeSum1E() = default;

    //! \brief The Bloch phase of an INTEGER cell offset \f$n\f$: \f$e^{ik\cdot R_n}\f$.  The k-CONVENTION stays
    //! entirely on the lattice/GPW side (the caller supplies this closure); the molecular side only ever asks
    //! "what is the phase of offset \f$n\f$" for the offsets it enumerates internally (a pre-built
    //! \c (Rs,phases) pair cannot work: the offsets are magnitude-screened per pair, never a fixed
    //! caller-visible set).  \f$\Gamma\f$ = the constant 1.
    using cellphase_t = std::function<dcmplx(const ivec3_t& n)>;

    //! \f$S_{ij}=\sum_R e^{ik\cdot R}\langle\chi_i|\chi_j(\cdot-R)\rangle\f$ (normalised), summed to
    //! \f$\varepsilon\f$ internally per shell pair.
    virtual chmat_t MakeOverlap(const cellphase_t& phase, const UnitCell& A) const = 0;

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
    //! \f$b_i=\sum_n \mathrm{phase}(n)\,\langle\chi_i|\,g(\cdot-C-R_n)\,\rangle\f$ -- the vector (\f$n\f$)
    //! analogue of the matrix \c MakeOverlap, with \a g standing in the \f$\chi_j\f$ slot; the same internal
    //! per-(\f$\chi_i\f$, \a g) magnitude-screened enumeration.  \f$\chi_i\f$ carries its usual
    //! normalisation; \a g is integrated RAW (its scale is in \c terms).
    virtual cvec_t MakeOverlap(const cellphase_t& phase, const UnitCell& A,
                               const GaussianFunction& g) const = 0;
    //! \brief The FINITE (home-term-only) overlap \f$b_i=\langle\chi_i|g\rangle\f$ -- the molecule/box limit
    //! (no lattice), used by the finite-mode KB assembly and any finite \f$\langle\chi|g\rangle\f$ consumer.
    virtual cvec_t MakeOverlap(const GaussianFunction& g) const = 0;
    //! \f$\langle p^2\rangle_{ij}=\sum_R e^{ik\cdot R}\langle\chi_i|-\nabla^2|\chi_j(\cdot-R)\rangle\f$ (NO 1/2 --
    //! the Hamiltonian applies it, matching \c Integrals_Kinetic).
    virtual chmat_t MakeKinetic(const cellphase_t& phase, const UnitCell& A) const = 0;
    //! \f$\sum_R e^{ik\cdot R}\langle\chi_i|\sum_c -Z_c/|r-R_c|\,|\chi_j(\cdot-R)\rangle\f$ over the atoms of \a cl.
    virtual chmat_t MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const = 0;

    //! \brief The lattice-summed matrix of a per-atom LOCAL Gaussian(-polynomial) operator placed at every atom
    //! of \a cl: \f$\langle\chi_i^k|\sum_a g_a|\chi_j^k\rangle=\sum_R e^{ik\cdot R}\sum_a\langle\chi_i|g_a|
    //! \chi_j(\cdot-R)\rangle\f$, the analytic 3-centre (\c Overlap3C) MATRIX sibling of the 2-centre
    //! \c MakeOverlap(g) VECTOR -- \a g stands in an operator slot BETWEEN the two orbitals (a multiplicative
    //! local potential), not the \f$\chi_j\f$ slot.  \a opForZ supplies the operator's \f$\{\alpha,\text{terms}\}\f$
    //! for a nuclear species \a Z (the impl centres it at each atom).  For a COMPACT operator (Gaussian tail,
    //! e.g. the short-range local pseudopotential) an image-atom operator is negligible (\f$e^{-|cell|^2/2r^2}\f$),
    //! so only home-cell atoms are placed -- like \c MakeNuclear, but here that is exact by construction.  The
    //! \f$\chi_i,\chi_j\f$ enumeration reuses the pair magnitude screen (an operator far from a screened-in pair
    //! contributes \f$\approx0\f$ via \c Overlap3C).  Hermitian; real at \f$\Gamma\f$.
    virtual chmat_t MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                                      const std::function<GaussianFunction(int Z)>& opForZ) const = 0;
    //! \brief The FINITE (home-term-only) \f$\langle\chi_i|\sum_a g_a|\chi_j\rangle\f$ -- the molecule/box limit
    //! (no lattice), the sibling of the finite \c MakeOverlap(g) used by the home-only GPW mode.
    virtual chmat_t MakeLocalGaussian(const Structure* cl,
                                      const std::function<GaussianFunction(int Z)>& opForZ) const = 0;

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

    //! \brief The STIFFNESS of the internal pair\f$\to\f$level assignment (the CP2K \c REL_CUTOFF safety): a
    //! pair \f$(i,j)\f$ demands a level with \f$e_{cut}\ge\f$ \c RelCutoffSafety()
    //! \f$\cdot\,e_{cut}^{ref}\,(\alpha_i+\alpha_j)/(2\alpha_{\max})\f$, so the stiffest possible requirement
    //! (the \f$\alpha_{\max}+\alpha_{\max}\f$ pair) is \c RelCutoffSafety() \f$\cdot\,e_{cut}^{ref}\f$.  A
    //! ladder is COMPLETE (every pair's requirement satisfiable) iff it contains a level at that cutoff --
    //! exposed as a scalar summary (like \c MaxExponent) so the ladder BUILDER can append the completion rung
    //! without duplicating the constant; the assignment itself stays internal (doc/GPWPlan.md 0b').
    virtual double RelCutoffSafety() const = 0;

    //! \brief ANALYTIC per-pair density collocation on a MULTI-GRID ladder (CP2K GPW / Quickstep): the grid
    //! density \f$\rho(r)=\sum_{ij}\sum_{R}\,\mathrm{Re}[D_{ij}e^{-ik\cdot R}]\,\chi_i(r)\chi_j(r-R)\f$ --
    //! the cell-periodic density of the Bloch orbital products \f$\chi_i^k\overline{\chi_j^k}\f$ -- with each
    //! (pair, cross-cell offset \f$R\f$) term evaluated ANALYTICALLY on its compact exp-tail box and
    //! MODULO-WRAPPED onto the grid.  The offsets are magnitude-screened (no hard \c Rcut, no Gibbs ringing);
    //! the wrap IS the image sum.  Density-matrix driven, so \f$\int\rho=\mathrm{Tr}(DS^k)\f$ (\f$S^k\f$ the
    //! screened-complete Bloch overlap).  MULTI-GRID (CP2K \c REL_CUTOFF): \a ecut_L[0] is the RESOLUTION
    //! REFERENCE (the charge-calibrated density grid); the remaining levels are conventionally coarser in
    //! descending order, and the ladder MAY append one FINER completion rung at
    //! \c RelCutoffSafety()\f$\cdot\f$\a ecut_L[0] (level selection is order-free; doc/GPWPlan.md 0b') --
    //! \a N_L[L] the level's grid divisions of cell \a A (raster \f$r=A(idx/N)\f$), \a ecut_L[L] its cutoff;
    //! each pair is collocated on the COARSEST level that still resolves its product exponent
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
    //! \a screenD: OPTIONAL density-magnitude screen -- when the caller supplies the density matrix whose
    //! field \a V_L is being integrated, each (pair, image) term is kept exactly when \c CollocateDensity of
    //! that density keeps it, so the collocate/integrate ADJOINT is exact on the shared truncated operator
    //! and the sweep only touches terms the density resolves (the CP2K eps/|coef| radii).  Density language
    //! only: this face already speaks \c chmat_t densities (\c CollocateDensity).
    virtual chmat_t IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                                       const std::vector<ivec3_t>& N_L,
                                       const std::vector<double>& ecut_L, double relCutoffScale=1.0,
                                       const chmat_t* screenD=nullptr) const = 0;
};

} //namespace
