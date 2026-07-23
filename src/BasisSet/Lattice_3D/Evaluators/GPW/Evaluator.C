// File: BasisSet/Lattice_3D/Evaluators/GPW/Evaluator.C  Gaussian-And-Plane-Waves orbital evaluator.
//
// The GPW sibling of PW_Evaluator (BasisSet/Lattice_3D/Evaluators/PW): a periodic ORBITAL evaluator whose
// orbitals are GAUSSIANS (compact, good for core/valence) instead of plane waves, standing on a lattice.  It
// satisfies the SAME isPW_1E_Evaluator concept, so the SAME EPW_Orbital1E_IBS<E> mixin builds the concrete
// GPW_IBS -- "GPW is a new evaluator, not a new IBS" (doc/MolecularPP_HarmonizationRound2.md section 2.5).
//
// The genuinely-new numerics of GPW's 1E tier are the PERIODIC Gaussian two-centre integrals -- the Bloch
// lattice sums  M_ij = Sum_R e^{ik.R} <chi_i | O | chi_j(.-R)>.  Rather than re-derive them, this evaluator
// DELEGATES them to the molecular Gaussian basis, which already owns the analytic McMurchie-Davidson kernels
// and GaussianRF::AtCenter (place a radial at an image): it holds a molecular orbital basis built over the
// cell's atoms and asks it the ONE high-level question Molecule::LatticeSum1E poses -- "sum these two-centre
// integrals over these lattice translations."  The Gaussian internals stay encapsulated on the molecular
// side; this evaluator only supplies the lattice geometry (the translation set from the cell) and complexifies
// the result.  No exponents or radials ever cross into qcLattice_BS.
//
// SCOPE (first increment): the GAMMA point (k=0).  The lattice sums are real (rsmat_t from the molecular side,
// widened to chmat_t here); general-k Bloch phases, the DFT tier (isPW_DFT_Evaluator: Hartree/XC by
// collocation) and a rigorous periodic nuclear attraction (Ewald) are later increments.
module;
#include <cstddef>
#include <functional>
#include <iosfwd>    // std::ostream (ReportGrids -- the grid diagnostic print)
#include <memory>
#include <string>
#include <vector>
export module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.BasisSet.Lattice_3D.Evaluators.PW;   // isPW_1E/DFT_Evaluator + PW_Grid_Evaluator (GPW's density grid)
import qchem.BasisSet.Molecule.LatticeSum1E;       // Molecule::LatticeSum1E (the periodic-1E capability we call)
import qchem.Pseudopotential.LocalPotential;       // LocalPotential_R  (the real-space local PP field V_loc(r))
import qchem.Pseudopotential.SeparablePotential;   // SeparablePotential_R (the KB projector radials beta_p(r))
import qchem.Mesh;                                 // qcMesh::MeshParams (the PP-quadrature integration mesh)
import qchem.BasisSet;                             // Real_BS (the molecular Gaussian basis we own)
import qchem.BasisSet.Orbital_1E_IBS;              // Real_OIBS (its orbital block: op()/Gradient/size)
export import qchem.BasisSet.Internal.GMap;        // G_ERI3 (the DFT 3-centre tensor, now with GPW weights)
import qchem.UnitCell;                             // UnitCell (the direct lattice; CellsInSphere/ToCartesian)
import qchem.Structure;                            // Structure (the nuclear-attraction centres)
import qchem.Types;                                // rvec3_t, cvec_t, cvec3vec_t, chmat_t, rmat_t

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The GPW orbital evaluator: Gaussian orbitals on a lattice at a k-point (Gamma this increment).  Owns
//! the molecular Gaussian basis built over the cell's atoms and reaches its periodic-1E capability
//! (\c Molecule::LatticeSum1E) by an abstract cross-cast; supplies the lattice translation set from the cell.
//! Satisfies \c isPW_1E_Evaluator, so the \c EPW_Orbital1E_IBS mixin drives the concrete \c GPW_IBS.
class GPW_Evaluator
{
public:
    //! \param mol   the molecular Gaussian orbital basis built over \a cell's atoms (kept alive here); its one
    //!              orbital block must realise \c Molecule::LatticeSum1E (a PG Gaussian basis does).
    //! \param cell  the direct lattice (source of the cell geometry + the reciprocal cell).
    //! \param densityEcut  the DENSITY/collocation grid selector (Hartree) -- GPW's ONLY grid cutoff (there is
    //!              no orbital/wavefunction cutoff; Gaussians are analytic).  THREE modes: \f$<0\f$ = AUTOMATIC
    //!              (the recommended default): the grid is floored to \a cutoffFactor\f$\cdot\alpha_{\max}\f$ from
    //!              the basis, so the caller need not know the Hartree value.  \f$=0\f$ = DFT tier OFF (1E-only).
    //!              \f$>0\f$ = EXPLICIT Hartree cutoff, honoured as given but with a \c cerr warning if it is below
    //!              \a cutoffFactor\f$\cdot\alpha_{\max}\f$ (under-resolves the density -> charge leaks off-grid).
    //! \param cutoffFactor  \f$C\f$ in the density-grid cutoff \f$E_{cut}=C\cdot\alpha_{\max}\f$ (default 2).  A
    //!              DENSITY-scale constant: the density is the PRODUCT of two orbitals, sharpest exponent
    //!              \f$2\alpha_{\max}\f$ — \f$C=2\f$ resolves the density at its own exponent.  HISTORY
    //!              (doc/GPWPlan 0.5(f)): the old default 8 was COMPENSATION for the ball-projected XC feed
    //!              (its Gibbs lobes went negative on sharp products and the \f$\rho>0\f$ guard collapsed
    //!              \f$E_{xc}\f$; only \f$8\alpha_{\max}\f$ ran clean).  The 0.5(f2) RAW collocated feed
    //!              (\f$\rho_{DM}\ge 0\f$ by construction) removed that failure mode, and the measured NaF
    //!              SR2 curve sits on a sub-2-mHa plateau from \f$C\approx 1.5\f$ up (Ecut=80/C=2: 1.2 mHa
    //!              from the C=8 answer, 0.1 mHa from CP2K).  CP2K's own operating point is \f$\sim 4\f$ in
    //!              these units via real-space collocation; our raw feed reaches lower because the finest
    //!              ladder level is sampled analytically, never truncated.
    //! \param kFrac fractional crystal momentum (fractional reciprocal coords; \f$\Gamma=0\f$).
    //! \param homeCellOnly  the FINITE-molecule MODE: no lattice images anywhere (1E matrices == the finite
    //!              molecule's; KB bra = the raw home orbital) -- the molecule-in-a-periodic-box configuration
    //!              the finite==lattice gates compare against.  A MODE, not a radius: in the periodic mode
    //!              (default) every lattice sum is an \f$\varepsilon\f$-CONVERGED SERIES enumerated internally
    //!              per shell pair -- THERE IS NO CUT in the R direction, and no radius parameter exists
    //!              (user pin, doc/GPWPlan.md).
    GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                  double densityEcut = 0.0, const rvec3_t& kFrac = rvec3_t(0,0,0),
                  bool homeCellOnly = false, double cutoffFactor = 2.0);
    //! Polymorphic (reached by the EPW_* mixin's Cast() cross-cast).  Releases this block's ladder-shaped
    //! collocation streams on the SHARED molecular evaluator (\c LatticeSum1E::ReleaseStreams) -- the streams
    //! are keyed by ladder shape, not by block, so without the release a finished stage's caches squat on the
    //! global stream budget for the process lifetime and starve any later grid (doc/GPWPlan.md 0.5(b): the
    //! grid-continuation fine stage at ~0% coverage).  Sibling k-blocks share the shape: the first dtor
    //! releases it, the rest no-op; a same-shape client that is still live simply rebuilds on demand.
    virtual ~GPW_Evaluator();

    // --- isPW_1E_Evaluator surface (exact signatures the concept demands) ---
    size_t     size()                          const {return itsN;}
    cvec_t     Eval(const rvec3_t& r)          const;   //!< \f$\sum_R e^{ik\cdot R}\chi_i(r-R)\f$ (real at \f$\Gamma\f$)
    cvec3vec_t EvalGradient(const rvec3_t& r)  const;   //!< \f$\sum_R e^{ik\cdot R}\nabla\chi_i(r-R)\f$
    chmat_t    OverlapMatrix()                 const;   //!< \f$\sum_R\langle i|j(\cdot-R)\rangle\f$
    chmat_t    KineticMatrix()                 const;   //!< \f$\sum_R\langle i|-\nabla^2|j(\cdot-R)\rangle\f$ (no 1/2)
    chmat_t    NuclearMatrix(const Structure* cl) const;//!< \f$\sum_R\langle i|\sum_c -Z_c/|r-R_c||j(\cdot-R)\rangle\f$

    // --- isPW_DFT_Evaluator surface: the DFT tier by COLLOCATION on the density grid (the genuinely-new GPW
    //     primitive).  All three route through the held PW_Grid_Evaluator (its {r}/{G}/FFT/Poisson). ---
    //! \brief Coulomb 3-centre tensor.  The per-column WEIGHT is the OVERLAP-metric (single-\f$r\f$) Fourier
    //! coefficient of the orbital product, \f$W_c(i,j)=\tfrac1\Omega\int\chi_i(r)\chi_j(r)\,e^{-iG_c\cdot r}\,d^3r\f$
    //! (\c ONE integration variable -- it is the collocated product FT, i.e. the density-side of the two-electron
    //! integral).  The SECOND electron coordinate and the \f$1/r_{12}\f$ live in the per-column diagonal Poisson
    //! kernel \f$4\pi/|G_c|^2\f$ (the \f$r_2\f$ integral \f$\int e^{iG_c\cdot r_2}/r_{12}\f$ in reciprocal space).
    //! So the full two-electron Coulomb is \c weight\f$\times\f$\c kernel, factorised through G-space -- never a
    //! 2-\f$r\f$ integral.  Density contracts \f$D\f$ against this via \c ContractG_ERI3.
    G_ERI3  Repulsion3CTensor() const;
    //! \brief Overlap 3-centre tensor: the same single-\f$r\f$ weight \f$W_c(i,j)\f$, EMPTY kernel -- the density's
    //! Fourier coefficient \f$\tilde\rho(G_c)=\sum_{ij}D_{ij}W_c(i,j)\f$ (no Poisson).
    G_ERI3  Overlap3CTensor() const;
    //! \brief The SAME two tables over an EXPLICIT fit grid \a grid (rather than the block's own \c DensityGrid).
    //! The tensor columns are \a grid's \f$\{G\}\f$ and the collocation multi-grid ladder derives from \a grid --
    //! so the fit-basis GRID/\f$\{G\}\f$ policy (set by \c CreateCD/VxcFitBasisSet) is HONOURED, not silently
    //! overridden by the block's own grid.  These are the seam \c GPW_IBS::MakeRepulsion3C/MakeOverlap3C call
    //! with the requested fit basis's grid; the no-arg overloads above are \c Repulsion3CTensor(itsFFT_R_G_Grids)
    //! convenience (tests + the block's own fit basis).  (doc/GPWPlan §0e -- "return the requested table".)
    G_ERI3  Repulsion3CTensor(std::shared_ptr<const PW_Grid_Evaluator> grid) const;
    G_ERI3  Overlap3CTensor  (std::shared_ptr<const PW_Grid_Evaluator> grid) const;
    //! \brief The potential->KS-matrix bridge (collocation's EXACT adjoint): \f$\langle\chi_i^k|V|\chi_j^k\rangle\f$
    //! by the ANALYTIC per-pair integrate-back (\c Molecule::LatticeSum1E::IntegratePotential) on the REL_CUTOFF
    //! multi-grid ladder -- \a Vtilde is restricted to each level's own \f$\{G\}\f$ (a SPECTRAL low-pass, no
    //! ringing), inverse-FFT'd to that level's grid, and each orbital pair gathers on the coarsest level that
    //! resolves its product exponent.  Only \f$V\f$ is ever sampled (weighted by the analytic Gaussians), so a
    //! diffuse pair on its matched coarse grid is accurate -- including against the SHARP local PP (the pair's
    //! own bandwidth bounds what it can sense of \f$V\f$).  Satisfies \c isPW_DFT_Evaluator; forwarded by
    //! \c EPW_Orbital_DFT_IBS to \c MakeOverlap.  (The CP2K analytic method, doc/GPWPlan.md \S0 -- replaces the
    //! sampled \c PhiOnGrid dense/patched/multigrid integrate-backs.)
    chmat_t OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    // --- Real-space external (pseudo)potential assembly: the GPW external term.  Unlike the plane-wave
    //     basis (G-space form factors, which Gaussians cannot supply) GPW quadratures the pseudopotential in
    //     REAL SPACE against its Gaussians on the cell's uniform integration mesh -- the SAME machinery the
    //     molecular PP_Local/PP_NonLocal terms use (qcMesh::WeightedOverlap / Overlap), so a Gaussian-in-a-box
    //     GPW run reproduces the finite molecular PP matrices.  At \f$\Gamma\f$ the matrices are real (widened
    //     to complex).  These realise Integrals_Pseudo<dcmplx> on GPW_IBS -> the whole Ham_PW_DFT drives GPW.
    //! \brief Local PP matrix, assembled in G-SPACE from the analytic form factor -- IDENTICALLY to the
    //! plane-wave path (\c PW_Evaluator::LocalPotentialMatrix), so GPW inherits the PW G=0 / alignment
    //! convention exactly.  \f$\langle\chi_i|V_{loc}|\chi_j\rangle=\langle\chi_i|\tilde V|\chi_j\rangle\f$ with
    //! \f$\tilde V(\Delta G)=\tfrac1\Omega\sum_a v_{loc}(Z_a,|\Delta G|^2)e^{-i\Delta G\cdot\tau_a}\f$,
    //! \f$\Delta G=0\f$ DROPPED (the ill-defined \f$-Z_{ion}/r\f$ cell-mean; carried by \c FormFactorG0
    //! alignment + Ewald background).  Assembled via \c OverlapMatrix (the collocation adjoint reconstructs
    //! \f$V_{loc}(r)\f$ on the density grid, band-limited to \c densityEcut, then quadratures it) -- so it is
    //! box-INDEPENDENT (unlike a real-space quadrature of the raw \f$V_{loc}\f$, whose Coulomb-tail cell-mean
    //! drifts with box size).  The KB nonlocal stays real-space (localized, no G=0 issue): see \c MakeSeparablePP.
    //! Which piece of the local PP a \c MakeLocalPP sweep assembles (the CP2K split, doc/GPWPlan.md 0e-PP):
    //! the FULL \f$V_{loc}\f$, its LONG (softened-Coulomb / deep-well) part, or its SHORT (compact
    //! poly-Gaussian) remainder.  All three route through the SAME sharp-field sweep, so \c Short+\c Long ==
    //! \c Full matrix-for-matrix (the sweep is linear in the form factor) -- the split relocates the LONG
    //! part's ENERGY into the Hartree term without changing the assembled matrices.  The sweep runs on the
    //! FULL ladder with the ABSOLUTE pair->level rule \f$e_{cut}\ge\kappa(\alpha_i+\alpha_j)\f$
    //! (\f$\kappa=30\f$ Ha, \f$e^{-\kappa/2}\f$ pair tails -- CP2K's \c REL_CUTOFF), which makes each
    //! piece STANDALONE-exact (no long/short cancellation partner needed; doc/GPWPlan.md 0e-PP step (a)).
    enum class LocalPart { Full, Long, Short };
    chmat_t MakeLocalPP    (const Structure* cl, const Pseudopotential::LocalPotential& loc,
                            LocalPart part=LocalPart::Full) const;
    //! \brief The LONG-range (softened-Coulomb) local-PP matrix \f$\langle i|V_{long}|j\rangle\f$, folded into
    //! the Hartree Poisson by \c PW_Hartree (doc/GPWPlan.md 0e-PP).  Convenience for \c MakeLocalPP(Long).
    chmat_t MakeLocalPPLong(const Structure* cl, const Pseudopotential::LocalPotential& loc) const;
    //! \brief The SHORT-range (compact poly \f$\times\f$ Gaussian) local-PP matrix \f$\langle i|V_{short}|j\rangle\f$.
    //! ANALYTIC when the model exposes its short part in closed Gaussian form (\c LocalPotential_Gaussian): a
    //! lattice-summed 3-centre \c Overlap3C (\c LatticeSum1E::MakeLocalGaussian), NO grid sweep -- the
    //! increment-2 cost win (doc/GPWPlan.md 0e-PP).  Falls back to the grid sweep (\c MakeLocalPP(Short)) for a
    //! model without the closed-Gaussian face.
    chmat_t MakeLocalPPShort(const Structure* cl, const Pseudopotential::LocalPotential& loc) const;
    //! \brief KB separable nonlocal matrix \f$\sum_{a,p,m}D_p|b\rangle\langle b|\f$ with the projection vector
    //! \f$b_i=\langle\chi_i|\beta_p(|r-R_a|)Y_{lm}\rangle\f$ (mesh quadrature).  Real symmetric at \f$\Gamma\f$.
    chmat_t MakeSeparablePP(const Structure* cl, const Pseudopotential::SeparablePotential_R& sep) const;

    //! \brief The density/collocation grid engine — ONE object carrying TWO layers (don't let the member
    //! name mislead): it IS-A \c PW_Evaluator = the Ecut BALL \f$\{G:\tfrac12|G|^2<E_{cut}\}\f$ (its
    //! \c size()/\c Gs(), i.e. the FIT-BASIS dimension \f$n_G\f$), and it HOLDS the FFT raster \f$N\f$
    //! (5-smooth-padded, alias-free) as that ball's \f$\{r\}\leftrightarrow\{G\}\f$ QUADRATURE engine.
    //! \c CreateCDFitBasisSet wraps this, so \f$\{G\}_\rho\f$ == the BALL; the raster is scaffolding
    //! (N is never a physics dial — the 5-smooth flip re-pinned nothing).
    const PW_Grid_Evaluator& DensityGrid() const {return *itsFFT_R_G_Grids;}

    //! \brief GRID DIAGNOSTIC (doc/GPWPlan §0e): the orbital-basis exponent line (\f$\alpha_{\min}/\alpha_{\max}\f$,
    //! \c cutoffFactor -- so the \f$C\cdot\alpha_{\max}\f$ grid policy is visible) followed by one line per STORED
    //! grid -- the FFT \f$\{r\}\leftrightarrow\{G\}\f$ density grid and every REL_CUTOFF ladder level (base
    //! sub-ladder = the local-PP integration grids; top completion rung flagged) -- each with FFT divisions
    //! \f$N\f$, \c Ecut, \f$n_G\f$ and the spectral extent \f$|G|_{\min}/|G|_{\max}\f$.  Printed at the start of
    //! every run (the basis-set ctor) so GPW's grids can be lined up against CP2K's \c &MGRID log output.
    void ReportGrids(std::ostream& os) const;
    //! One diagnostic line for any stored grid (the per-grid piece of \c ReportGrids; also used by the
    //! \c CreateCD/VxcFitBasisSet factories to print the fit grid they actually wrap).
    static std::ostream& ReportGrid(std::ostream& os, const std::string& tag, const PW_Grid_Evaluator& g);

    //! Cache-key fragment: the molecular basis's ID + \f$k\f$ + translation count + the density-grid cutoff
    //! (the collocation tensor depends on the grid, so the framework cache key must pin it).
    std::string IDFragment() const;

private:
    std::shared_ptr<const BasisSet::Real_BS> itsMol;   //!< owns the molecular Gaussian basis (lifetime)
    const BasisSet::Real_OIBS*          itsOrb = nullptr; //!< its single orbital block (op()/Gradient/size)
    const Molecule::LatticeSum1E*       itsLat = nullptr; //!< the same block's periodic-1E capability (cross-cast)
    // INTERNAL Bloch-orbital image set for Eval/EvalGradient + the mesh-path KB quadrature (the ONLY places
    // that still walk an explicit image list -- the 1E matrices and the analytic KB/collocation enumerate
    // inside the molecular seam).  DERIVED from the eps-screen (single-orbital reach + cell span; home-only
    // mode = the origin alone) -- an implementation detail of the series' convergence, never a parameter.
    std::vector<rvec3_t>                itsRc;            //!< Eval/mesh-KB translations (incl. origin)
    cvec_t                              itsPhaseC;        //!< matching Bloch phases e^{ik.R} (origin = 1)
    bool                                itsHomeOnly=false; //!< the finite-molecule MODE (no images anywhere)
    rvec3_t                             itsk;             //!< fractional crystal momentum k
    double  itsMaxReach=0.0;   //!< max orbital reach sqrt(-ln eps/alpha_min) (Eval per-image screen)
    rvec3_t itsCellCtr;        //!< cell centre + bounding radius: an image R contributes at r only if
    double  itsCellRad=0.0;    //!< |(r-R)-ctr| <= rad+maxReach (every orbital centre sits inside the cell)
    UnitCell                            itsCell;          //!< the DIRECT cell (stored, NOT reconstructed by a
                                                          //!< double-reciprocal round trip -- that reconstruction
                                                          //!< can re-orient a non-cubic (FCC) primitive cell and
                                                          //!< silently shift every collocation box)
    size_t                              itsN   = 0;       //!< number of Gaussian orbitals
    double  itsCutoffFactor=2.0;   //!< the density-grid floor constant C (ctor param) -- the ENERGY calibration
                                   //!< C*RelCutoffSafety()*alpha_max gates the top completion rung (EnsureLevels)
    std::shared_ptr<const PW_Grid_Evaluator> itsFFT_R_G_Grids;     //!< the density/collocation grid (null if DFT tier off)
    // NO hand-rolled tensor cache: the collocation tensor is a stateless build; the FRAMEWORK caches it
    // (BasisSet::Band_FT_IBS::Repulsion3C/Overlap3C via theCache<dcmplx>(), keyed by BasisSetID -- see IDFragment).
    //! \brief Last-density collocation memo, SHARED by the two \c MakeCollocator closures (Coulomb + overlap
    //! tensors): each SCF iteration collocates the SAME \f$D\f$ twice -- once for \f$\rho\f$ on the grid
    //! (XC/charge) and once for the Coulomb apply -- so the second call replays the level densities for free
    //! (exact-equality keyed on \f$D\f$; the phase/cell/ladder are fixed per evaluator, i.e. per k-block).
    //! The G-space post-processing (FFT + nested combine +/- Coulomb kernel) stays per-call.
    //! \c ecut pins the LADDER the cached \c rho was collocated on: the memo dedups the two same-grid closures
    //! (Coulomb + overlap) collocating the same \f$D\f$ per iteration, but the SAME \f$D\f$ can be collocated on
    //! a DIFFERENT grid (a coarse block asked for a fine fit grid -- the grid-continuation seed), where replaying
    //! the wrong-ladder \c rho is a dimension mismatch (segfault).  So \c sameD requires \c ecut to match too.
    struct CollocMemo { bool valid=false; chmat_t D; std::vector<rvec_t> rho; std::vector<double> ecut; };
    mutable std::shared_ptr<CollocMemo> itsCollocMemo;    //!< shared so both framework-cached closures see it
    //! The Bloch phase of an integer cell offset \f$n\f$: \f$e^{2\pi i\,k_{frac}\cdot n}\f$ -- the closure the
    //! analytic kernels call back for each screened cross-cell pair offset (the k-CONVENTION stays here,
    //! lattice-side; the molecular basis never sees \f$k\f$).
    Molecule::LatticeSum1E::cellphase_t CellPhase() const;
    //! The matrix-free density->\f$\tilde\rho\f$ (\a coulomb=false) or ->\f$V_H\f$ (\a coulomb=true,
    //! \f$4\pi/G^2\f$ folded in) map, as a closure: ANALYTIC per-pair collocation on the multi-grid ladder
    //! (\c LatticeSum1E::CollocateDensity), one FFT per level, \f$\tilde\rho\f$ combined NESTED in G-space --
    //! the \c G_ERI3::apply realization.
    //! The collocator over an EXPLICIT fit grid: the ladder derives from \a grid (the tensor's requested
    //! \f$\{G\}\f$), NOT the block's own \c itsFFT_R_G_Grids -- the closure captures its own ladder built from \a grid.
    std::function<ΔG_Map(const chmat_t&)> MakeCollocator(bool coulomb, std::shared_ptr<const PW_Grid_Evaluator> grid) const;
    //! The BACKWARD adjoint of \c MakeCollocator: a self-contained field->matrix integrate-back over \a grid's
    //! ladder -- the \c G_ERI3::applyAdjoint the Overlap3C/Repulsion3C(grid) tensors carry (doc/GPWPlan §0e step2).
    std::function<chmat_t(const std::function<dcmplx(const ivec3_t&)>&)>
        MakeIntegrator(std::shared_ptr<const PW_Grid_Evaluator> grid) const;
    //! \brief The RAW-RASTER pair (doc/GPWPlan 0.5(f2)) -- \c G_ERI3::applyRaw / \c applyRawAdjoint.
    //! \c MakeRawCollocator: \f$D\to\rho_{DM}(r)\f$ on \a grid's raster -- the finest ladder level kept RAW
    //! (its analytic samples are pointwise \f$\ge 0\f$ to screening-\f$\varepsilon\f$ for PSD \f$D\f$), every
    //! other level transferred SPECTRALLY (zero-pad up / band-drop down, no ball restriction anywhere).
    //! \c MakeRawIntegrator: the exact transpose -- \f$v(r)\f$ band-truncates to each level's box, inverse-FFTs,
    //! and gathers through the SAME analytic \c IntegratePotential, so \f$H_{xc}=\partial E_{xc}/\partial D\f$
    //! exactly (the FFT normalisations cancel against the \f$\Omega/N_{pts}\f$ quadrature weights -- no factors).
    //! Both share \c CollocMemo with the ball closures (one collocation per (D, ladder) per iteration).
    std::function<rvec_t(const chmat_t&)> MakeRawCollocator(std::shared_ptr<const PW_Grid_Evaluator> grid) const;
    std::function<chmat_t(const rvec_t&)> MakeRawIntegrator(std::shared_ptr<const PW_Grid_Evaluator> grid) const;
    qcMesh::MeshParams PPMeshParams() const;  //!< the PP-quadrature integration mesh params (uniform, eCut=densityEcut)

    // The REL_CUTOFF multi-grid level ladder: the fine density grid + coarser grids (a factor 4 in Ecut each)
    // down to the level resolving the most-diffuse pair product (~Ecut*alpha_min/alpha_max), PLUS the TOP
    // COMPLETION RUNG at RelCutoffSafety()*Ecut appended LAST (doc/GPWPlan.md 0b': the sharpest pairs'
    // requirement exceeds the reference grid by construction; the rung makes every smooth-path assignment
    // satisfiable).  ecut_L[0] stays the RESOLUTION REFERENCE (the density grid) -- level selection on the
    // molecular side is order-free.  Built once (geometry-fixed) by EnsureLevels.  The density collocation
    // (MakeCollocator), the integrate-back (OverlapMatrix) and the local-PP sweep (MakeLocalPP, absolute
    // kappa rule) all run per-pair on the FULL ladder; itsNBaseLevels only labels the top rung (diagnostics).
    void EnsureLevels() const;   //!< builds/caches the block's OWN ladder (itsFFT_R_G_Grids) for OverlapMatrix + MakeLocalPP
    //! Build the REL_CUTOFF ladder from an ARBITRARY fine grid \a grid (level [0]==\a grid) into the output
    //! vectors -- the grid-parameterized core of \c EnsureLevels, so the collocator can build a ladder from the
    //! REQUESTED fit grid while the block's own paths keep the cached \c itsLevels.
    void BuildLevels(std::shared_ptr<const PW_Grid_Evaluator> grid,
                     std::vector<std::shared_ptr<const PW_Grid_Evaluator>>& levels,
                     std::vector<ivec3_t>& levelN, std::vector<double>& levelEcut, size_t& nBaseLevels) const;
    mutable std::vector<std::shared_ptr<const PW_Grid_Evaluator>> itsLevels;   //!< [0]==itsFFT_R_G_Grids (reference); coarser; top rung LAST
    mutable std::vector<ivec3_t>    itsLevelN;     //!< each level's FFT grid divisions
    mutable std::vector<double>     itsLevelEcut;  //!< each level's cutoff (reference first)
    mutable size_t                  itsNBaseLevels=0; //!< levels before the top rung (the local-PP sub-ladder)
};

static_assert(isPW_1E_Evaluator <GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_1E_Evaluator");
static_assert(isPW_DFT_Evaluator<GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_DFT_Evaluator");

} //namespace
