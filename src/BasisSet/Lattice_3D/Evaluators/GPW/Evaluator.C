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
    //! \param Rcut  radius (a.u.) of the OVERLAP/1E lattice-translation sphere; \f$\le 0\f$ means the home cell
    //!              only (\f$R=0\f$: reproduces the finite molecule exactly).  The analytic 2-centre lattice
    //!              sums (overlap/kinetic/nuclear) need a generous \a Rcut to converge to a positive-definite
    //!              overlap (the truncated single sum is indefinite at small \a Rcut).
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
    //! \brief The potential->KS-matrix bridge (collocation's adjoint): \f$\langle\chi_i|V|\chi_j\rangle=\int\chi_i
    //! V\chi_j\f$ with \f$V(r)\f$ the inverse-FFT of \a Vtilde over the density grid -- grid-integrate, not the PW
    //! Fourier lookup.  Satisfies \c isPW_DFT_Evaluator; forwarded by \c EPW_Orbital_DFT_IBS to \c MakeOverlap.
    chmat_t OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;
    //! The dense (fine-grid) integrate-back \f$w\,\Phi^H(V\!\odot\!\Phi)\f$ via OpenBLAS zgemm -- the default
    //! path, and the one the SHARP static local PP always uses (\c MakeLocalPP).  \c OverlapMatrix dispatches here
    //! unless the multi-grid path is enabled (\c UseMultiGrid), which routes only the smooth dynamic V_H+V_xc.
    chmat_t DenseOverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    //! \brief The PATCHED integrate-back: same result as \c OverlapMatrix(Vtilde) but delegated to the
    //! molecular-side \c Molecule::LatticeSum1E::MakePotentialMatrix -- per-orbital Gaussian-support patches
    //! (the localized part), contracting each pair only on its support overlap.  The Gaussian primitives stay
    //! encapsulated (only \f$V(r)\f$ + the grid cross the seam).  Single-grid scaffold for the multi-grid
    //! rewrite (doc/GPWPlan.md \S0); OPT-IN (the default \c OverlapMatrix stays the byte-identical dense GEMM),
    //! validated bit-consistent in \c GPW_UT.  Will become the default once the grid levels make it win.
    chmat_t PatchedOverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    //! \brief The MULTI-GRID integrate-back (doc/GPWPlan.md \S0 Increment 2 -- the NaF gap-closer): builds a
    //! ladder of density grids from the fine \c densityEcut down to \f$\propto\alpha_{\min}\f$ (each a factor 4
    //! in \f$E_{cut}\f$), restricts \a Vtilde to each level, and delegates to
    //! \c Molecule::LatticeSum1E::MakePotentialMatrixMG, which contracts each orbital pair on the COARSEST
    //! level resolving its product exponent.  Diffuse pairs then live on small coarse grids instead of the fine
    //! grid dictated by the tightest primitive.  APPROXIMATE vs the single fine grid (converges as the ladder
    //! refines); OPT-IN.  Reduces to \c PatchedOverlapMatrix when the basis spans a single exponent decade.
    chmat_t MultiGridOverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    //! Enable the multi-grid path for the DYNAMIC (per-iteration) integrate-back (\c OverlapMatrix(Vtilde)).
    //! The STATIC local PP stays dense (\c MakeLocalPP). Default OFF (dense) -- committed anchors byte-identical.
    void UseMultiGrid(bool on=true) const {itsUseMG=on;}

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
    size_t                              itsN   = 0;       //!< number of Gaussian orbitals
    std::shared_ptr<const PW_Grid_Evaluator> itsGrid;     //!< the density/collocation grid (null if DFT tier off)
    // NO hand-rolled tensor cache: the collocation tensor is a stateless build; the FRAMEWORK caches it
    // (BasisSet::Band_FT_IBS::Repulsion3C/Overlap3C via theCache<dcmplx>(), keyed by BasisSetID -- see IDFragment).
    mutable mat_t<dcmplx> itsPhiOnGrid;  //!< cached \f$\chi_i^k(r_g)\f$ (geometry-fixed; built once, reused every SCF iteration)
    //! Reused SCRATCH buffer (NOT a cache -- its contents \f$V\odot\Phi\f$ change every SCF iteration).  Sized
    //! once (Npts x n, ~130 MB for NaF); OverlapMatrix refills + reuses it so the zgemm operand is not
    //! re-allocated on every call (~120 calls/run -> otherwise ~16 GB of alloc/free churn + page faults).
    mutable mat_t<dcmplx> itsVPhiBuf;
    //! \f$\chi_i^k(r_g)\f$ on the density grid (Npts x n; real at \f$\Gamma\f$).  CACHED: it is a pure function
    //! of geometry (grid + orbitals + image set), identical across SCF iterations, and the per-iteration
    //! integrate-back \c OverlapMatrix reuses it -- recomputing the Bloch sum at every grid point per iteration
    //! was the dominant cost at a large image \c Rcut (the density-independent \f$O(N_{pts}\cdot|R|)\f$ sum).
    const mat_t<dcmplx>& PhiOnGrid() const;
    G_ERI3  BuildWeights() const;  //!< the collocation weight tensor (columns=grid {G}, weights, NO kernel)
    qcMesh::MeshParams PPMeshParams() const;  //!< the PP-quadrature integration mesh params (uniform, eCut=densityEcut)

    // MULTI-GRID (Increment 2) level ladder: the fine density grid + coarser grids (a factor 4 in Ecut each)
    // down to ~cutoffFactor*alpha_min, built once (geometry-fixed) by EnsureLevels and cached with their points /
    // cutoffs / quadrature weights.  MultiGridOverlapMatrix restricts Vtilde to each level and delegates the
    // per-pair contraction to the molecular side.  Empty unless the multi-grid path is used.
    void EnsureLevels() const;
    mutable std::vector<std::shared_ptr<const PW_Grid_Evaluator>> itsLevels;   //!< finest first; [0]==itsGrid
    mutable std::vector<rvec3vec_t> itsLevelPts;   //!< each level's grid points (geometry-fixed)
    mutable std::vector<double>     itsLevelEcut;  //!< each level's cutoff (descending)
    mutable std::vector<double>     itsLevelW;     //!< each level's quadrature weight Omega/Npts(L)
    mutable bool itsUseMG=false;     //!< opt-in: route the dynamic V_H+V_xc integrate-back via multi-grid (UseMultiGrid)
};

static_assert(isPW_1E_Evaluator <GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_1E_Evaluator");
static_assert(isPW_DFT_Evaluator<GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_DFT_Evaluator");

} //namespace
