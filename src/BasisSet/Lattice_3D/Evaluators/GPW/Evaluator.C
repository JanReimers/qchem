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
    //! \param cell  the direct lattice (source of the real-space translation set \f$\{R\}\f$ + the reciprocal cell).
    //! \param densityEcut  the DENSITY/collocation grid selector (Hartree) -- GPW's ONLY grid cutoff (there is
    //!              no orbital/wavefunction cutoff; Gaussians are analytic).  THREE modes: \f$<0\f$ = AUTOMATIC
    //!              (the recommended default): the grid is floored to \a cutoffFactor\f$\cdot\alpha_{\max}\f$ from
    //!              the basis, so the caller need not know the Hartree value.  \f$=0\f$ = DFT tier OFF (1E-only).
    //!              \f$>0\f$ = EXPLICIT Hartree cutoff, honoured as given but with a \c cerr warning if it is below
    //!              \a cutoffFactor\f$\cdot\alpha_{\max}\f$ (under-resolves the density -> charge leaks off-grid).
    //! \param cutoffFactor  \f$C\f$ in the density-grid floor \f$C\cdot\alpha_{\max}\f$ (default 4, the calibrated
    //!              minimum; pass \f$C\ge4\f$ for a finer grid).  This is a DENSITY-scale constant: the density is
    //!              the product of two orbitals (exponent \f$2\alpha_{\max}\f$), so \f$C\f$ already folds in the
    //!              \f$\times2\f$ over a single-orbital cutoff.  See doc/GPWPlan.md \S0.
    //! \param kFrac fractional crystal momentum (fractional reciprocal coords; \f$\Gamma=0\f$).
    //! \param Rcut  the OVERLAP/1E lattice-translation sphere (a.u.).  THREE modes (mirroring \a densityEcut):
    //!              \f$<0\f$ = AUTOMATIC (recommended) -- the radius is derived from the basis
    //!              (\f$2\sqrt{-\ln\varepsilon/\alpha_{\min}}\f$ + the cell span), i.e. everything the
    //!              magnitude screen can keep is enumerated and the screen prunes it sparse; \f$=0\f$ = home
    //!              cell only (the FINITE-molecule mode); \f$>0\f$ = explicit radius (legacy).
    //! \param collRcut  radius (a.u.) of the COLLOCATION sphere -- the reach of the Bloch orbital sum
    //!              \f$\sum_R e^{ik\cdot R}\chi(r-R)\f$ the density grid samples.  The density is local
    //!              (on-site + nearest images), so this is MUCH smaller than \a Rcut -- decoupling the two
    //!              image sets is what makes multi-k bulk affordable (the collocation re-sums images at every
    //!              grid point, per SCF iteration -- the dominant cost).  \f$\le 0\f$ reuses \a Rcut's set
    //!              (backward-compatible: \f$\Gamma\f$/finite runs collocate on the same home-cell set).
    GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                  double densityEcut = 0.0, const rvec3_t& kFrac = rvec3_t(0,0,0), double Rcut = 0.0,
                  double collRcut = 0.0, double cutoffFactor = 4.0);
    virtual ~GPW_Evaluator() = default;   // polymorphic: reached by the EPW_* mixin's Cast() cross-cast

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
    chmat_t MakeLocalPP    (const Structure* cl, const Pseudopotential::LocalPotential& loc) const;
    //! \brief KB separable nonlocal matrix \f$\sum_{a,p,m}D_p|b\rangle\langle b|\f$ with the projection vector
    //! \f$b_i=\langle\chi_i|\beta_p(|r-R_a|)Y_{lm}\rangle\f$ (mesh quadrature).  Real symmetric at \f$\Gamma\f$.
    chmat_t MakeSeparablePP(const Structure* cl, const Pseudopotential::SeparablePotential_R& sep) const;

    //! The density/collocation grid engine (the fit basis is built over it, so \f$\tilde\rho\f$'s \f$\{G\}\f$ matches).
    const PW_Grid_Evaluator& DensityGrid() const {return *itsGrid;}

    //! Cache-key fragment: the molecular basis's ID + \f$k\f$ + translation count + the density-grid cutoff
    //! (the collocation tensor depends on the grid, so the framework cache key must pin it).
    std::string IDFragment() const;

private:
    std::shared_ptr<const BasisSet::Real_BS> itsMol;   //!< owns the molecular Gaussian basis (lifetime)
    const BasisSet::Real_OIBS*          itsOrb = nullptr; //!< its single orbital block (op()/Gradient/size)
    const Molecule::LatticeSum1E*       itsLat = nullptr; //!< the same block's periodic-1E capability (cross-cast)
    // The lattice translations {R} and their Bloch phases {e^{ik.R}} are a weighted point set, built + kept
    // together (future: ONE qcMesh cMesh = Mesh<dcmplx>).  Both always include the origin (phase 1) first.
    // TWO decoupled sets: the OVERLAP/1E set (large Rcut -> PSD overlap) and the COLLOCATION set (small
    // collRcut -> the local density's orbital reach; == the overlap set when collRcut<=0).
    std::vector<rvec3_t>                itsR;             //!< overlap/1E translations (incl. origin)
    cvec_t                              itsPhase;         //!< matching Bloch phases e^{ik.R} (origin = 1)
    std::vector<rvec3_t>                itsRc;            //!< collocation translations (orbital reach; incl. origin)
    cvec_t                              itsPhaseC;        //!< matching collocation Bloch phases (origin = 1)
    rvec3_t                             itsk;             //!< fractional crystal momentum k
    double  itsMaxReach=0.0;   //!< max orbital reach sqrt(-ln eps/alpha_min) (Eval per-image screen)
    rvec3_t itsCellCtr;        //!< cell centre + bounding radius: an image R contributes at r only if
    double  itsCellRad=0.0;    //!< |(r-R)-ctr| <= rad+maxReach (every orbital centre sits inside the cell)
    UnitCell                            itsCell;          //!< the DIRECT cell (stored, NOT reconstructed by a
                                                          //!< double-reciprocal round trip -- that reconstruction
                                                          //!< can re-orient a non-cubic (FCC) primitive cell and
                                                          //!< silently shift every collocation box)
    size_t                              itsN   = 0;       //!< number of Gaussian orbitals
    double  itsCutoffFactor=4.0;   //!< the density-grid floor constant C (ctor param) -- the ENERGY calibration
                                   //!< C*RelCutoffSafety()*alpha_max gates the top completion rung (EnsureLevels)
    std::shared_ptr<const PW_Grid_Evaluator> itsGrid;     //!< the density/collocation grid (null if DFT tier off)
    // NO hand-rolled tensor cache: the collocation tensor is a stateless build; the FRAMEWORK caches it
    // (BasisSet::Band_FT_IBS::Repulsion3C/Overlap3C via theCache<dcmplx>(), keyed by BasisSetID -- see IDFragment).
    //! \brief Last-density collocation memo, SHARED by the two \c MakeCollocator closures (Coulomb + overlap
    //! tensors): each SCF iteration collocates the SAME \f$D\f$ twice -- once for \f$\rho\f$ on the grid
    //! (XC/charge) and once for the Coulomb apply -- so the second call replays the level densities for free
    //! (exact-equality keyed on \f$D\f$; the phase/cell/ladder are fixed per evaluator, i.e. per k-block).
    //! The G-space post-processing (FFT + nested combine +/- Coulomb kernel) stays per-call.
    struct CollocMemo { bool valid=false; chmat_t D; std::vector<rvec_t> rho; };
    mutable std::shared_ptr<CollocMemo> itsCollocMemo;    //!< shared so both framework-cached closures see it
    //! The Bloch phase of an integer cell offset \f$n\f$: \f$e^{2\pi i\,k_{frac}\cdot n}\f$ -- the closure the
    //! analytic kernels call back for each screened cross-cell pair offset (the k-CONVENTION stays here,
    //! lattice-side; the molecular basis never sees \f$k\f$).
    Molecule::LatticeSum1E::cellphase_t CellPhase() const;
    //! The matrix-free density->\f$\tilde\rho\f$ (\a coulomb=false) or ->\f$V_H\f$ (\a coulomb=true,
    //! \f$4\pi/G^2\f$ folded in) map, as a closure: ANALYTIC per-pair collocation on the multi-grid ladder
    //! (\c LatticeSum1E::CollocateDensity), one FFT per level, \f$\tilde\rho\f$ combined NESTED in G-space --
    //! the \c G_ERI3::apply realization.
    std::function<ΔG_Map(const chmat_t&)> MakeCollocator(bool coulomb) const;
    qcMesh::MeshParams PPMeshParams() const;  //!< the PP-quadrature integration mesh params (uniform, eCut=densityEcut)

    // The REL_CUTOFF multi-grid level ladder: the fine density grid + coarser grids (a factor 4 in Ecut each)
    // down to the level resolving the most-diffuse pair product (~Ecut*alpha_min/alpha_max), PLUS the TOP
    // COMPLETION RUNG at RelCutoffSafety()*Ecut appended LAST (doc/GPWPlan.md 0b': the sharpest pairs'
    // requirement exceeds the reference grid by construction; the rung makes every smooth-path assignment
    // satisfiable).  ecut_L[0] stays the RESOLUTION REFERENCE (the density grid) -- level selection on the
    // molecular side is order-free.  Built once (geometry-fixed) by EnsureLevels.  The density collocation
    // (MakeCollocator) and the integrate-back (OverlapMatrix) run per-pair on the FULL ladder; the SHARP-field
    // local PP (MakeLocalPP, relCutoffScale=6) uses the BASE sub-ladder only (itsNBaseLevels) -- its stiffened
    // requirement would flood the top rung with mid pairs whose boxes are huge on the doubled grid.
    void EnsureLevels() const;
    mutable std::vector<std::shared_ptr<const PW_Grid_Evaluator>> itsLevels;   //!< [0]==itsGrid (reference); coarser; top rung LAST
    mutable std::vector<ivec3_t>    itsLevelN;     //!< each level's FFT grid divisions
    mutable std::vector<double>     itsLevelEcut;  //!< each level's cutoff (reference first)
    mutable size_t                  itsNBaseLevels=0; //!< levels before the top rung (the local-PP sub-ladder)
};

static_assert(isPW_1E_Evaluator <GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_1E_Evaluator");
static_assert(isPW_DFT_Evaluator<GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_DFT_Evaluator");

} //namespace
