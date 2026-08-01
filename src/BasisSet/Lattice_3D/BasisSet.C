// File: BasisSet/Lattice_3D/BasisSet.C  Factory + container for 3D-periodic (plane-wave) basis sets.
//
// The public entry point for building a crystal basis set: hand it a Lattice_3D (the cell + its
// Brillouin-zone grid) and a cutoff, get back an abstract tBasisSet<dcmplx>.  The concrete container
// (PW_BasisSet) and the per-k PlaneWave_IBS list it owns are an implementation detail (Imp/BasisSet.C);
// callers are forced through the polymorphic BasisSet interface, exactly as for molecular bases.
module;
#include <memory>
#include <vector>
export module qchem.BasisSet.Lattice_3D.BasisSet;
export import qchem.BasisSet;                          // Complex_BS (= tBasisSet<dcmplx>)
export import qchem.Lattice_3D;                        // Lattice_3D (the crystal structure + BZ grid)
export import qchem.BasisSet.Lattice_3D.PlaneWave_IBS; // PlaneWave_IBS + LocalPotential/SeparablePotential
export import qchem.BasisSet.Lattice_3D.GPW_IBS;       // GPW_IBS + CellImages (the GPWFactory mode argument)
export import qchem.BasisSet.Lattice_3D.Evaluators.PW; // RasterPolicy (a PUBLIC factory knob since 0.5(a))
import qchem.BasisSet.Internal.BasisSetImp;            // BasisSetImp<dcmplx> (the PW_BasisSet base; NOT re-exported)
import qchem.Types;                                    // dcmplx

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief Which 3D-periodic basis to build.  PW = plane waves (lineage A); APW/LAPW (lineage B) follow.
enum class Type { PW };

//! \brief Build a 3D-periodic basis set for \a lat at energy cutoff \a Ecut.  Returns an abstract
//! tBasisSet<dcmplx> on the heap (caller owns), so callers must use the polymorphic interface.
//! \note The pseudopotential model is NOT configured here: it lives on the Vpseudo Hamiltonian term
//!   (the pseudo-wall), which calls the basis's MakeLocalPotential/MakeSeparablePotential assembly.
//! \note Single-k for now: the returned basis holds ONE Bloch block at \f$\Gamma\f$.  Phase 2
//!   generalises this to the full BZ k-list (intended to be the only k-loop in the framework).
Complex_BS* Factory(Type type, const ::qchem::Lattice_3D& lat, double Ecut);

//! \brief Build a GPW basis (periodic Gaussians on the lattice) over \a lat from a molecular Gaussian basis
//! \a mol built on the cell's atoms, at density-collocation cutoff \a densityEcut.
//! Returns an abstract tBasisSet<dcmplx> (caller owns).  Unlike the PW \c Factory this needs a Gaussian
//! orbital basis (GPW = Gaussian orbitals); the pseudopotential still lives on the Hamiltonian term, which
//! reaches GPW's real-space \c Integrals_Pseudo<dcmplx> assembly -- so the same \c Ham_PW_DFT drives it.
//! THERE IS NO CUT (doc/GPWPlan.md pin): every lattice sum is an eps-converged series enumerated inside the
//! molecular seam -- no radius parameter exists on this surface.
//! \param kShift  fractional Monkhorst-Pack offset of the k-mesh (\f$0\f$ = Γ-centred; \f$½\f$ = the classic MP
//!                offset, i.e. CP2K's default for even grids -- \f$k=\pm¼\f$ at \f$N=2\f$).
//! \param densityEcut  \f$<0\f$ = AUTOMATIC density grid \a cutoffFactor\f$\cdot\alpha_{\max}\f$ (recommended);
//!        \f$=0\f$ = 1E-only; \f$>0\f$ = explicit Hartree cutoff (\c cerr warning if under-resolved).
//! \param images  the lattice-image MODE (\c CellImages::Periodic default; \c HomeCellOnly = the box gates).
//! \param cutoffFactor  \f$C\ge4\f$ in the density-grid floor \f$C\cdot\alpha_{\max}\f$ (default 4).
Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, rvec3_t kShift={0,0,0},
                       CellImages images=CellImages::Periodic, double cutoffFactor=2.0);

//! \brief The GPW factory knobs as ONE readable struct (doc/GPWPlan1.md item 1) -- designated
//! initializers make test recipes self-documenting: GPWFactory(lat, mol, {.densityEcut=40.0}).
//! STANDARD (the one users touch): \c densityEcut.  ADVANCED (sensible defaults, rarely changed):
//! everything else -- the expert-system goal is that the guards make hand-tuning unnecessary.
struct GPWParams
{
    // --- Standard ---
    double       densityEcut = -1.0;    //!< <0 AUTO floor cutoffFactor*alpha_max (recommended); =0 1E-only; >0 explicit Ha
    // --- Advanced ---
    double       cutoffFactor= 2.0;     //!< C in the auto floor C*alpha_max (2 = the density's own product exponent)
    RasterPolicy raster = RasterPolicy::BallOnly;   //!< FFT-raster policy.  DEFAULT BallOnly (user 2026-07-23,
                                                    //!< after the 0.5(a) A/B: ~1 mHa at/above the C=2 floor, ~8x
                                                    //!< fewer raster points; AliasFree = the exact-quadrature
                                                    //!< option -- kernel gates and sub-floor grids want it)
    CellImages   images = CellImages::Periodic;     //!< lattice-image MODE (HomeCellOnly = the finite-box gates)
    rvec3_t      kShift = rvec3_t(0,0,0);           //!< fractional MP offset of the k-mesh (0 = Gamma-centred)
    double       ladderFactor = 4.0;                //!< REL_CUTOFF multigrid PROGRESSION factor (CP2K's, default 3;
                                                    //!< ours 4): the per-step Ecut ratio of the coarse-level ladder.
                                                    //!< The DEPTH is automatic (~log_f(alpha_max/alpha_min) levels,
                                                    //!< cell-capped) -- this only tunes the GRADATION.  Smaller =
                                                    //!< more, finer-spaced levels (rarely needed; the auto ladder
                                                    //!< already routes each diffuse pair to a matched coarse grid).
    bool         reduceBZ = false;                  //!< IBZ/k-star: fold the MP mesh to its irreducible wedge under
                                                    //!< the crystal point group (τ=0 symmorphic ops only -- the safe
                                                    //!< guard; non-symmorphic crystals fold less but stay correct).
                                                    //!< Only the linear parts are used, so the density MUST be
                                                    //!< symmetrized over the star for this to be exact (doc/GPWPlan1
                                                    //!< item 3 IBZ).  DEFAULT off (full mesh, no symmetrization).
                                                    //!< This IS the §3 SymmetryPolicy assert-switch
                                                    //!< (doc/SymmetryUpgradePlan.md): true = the caller ASSERTS
                                                    //!< imposition of the full detected group (release-check audits
                                                    //!< SSB); false = a FREE run (nothing imposed; the detected group
                                                    //!< still feeds the order-parameter diagnostic).
    RasterFields rasterFields = RasterFields::HartreeXC;  //!< field-sharpness ROUTING: which terms the raster serves.
                                                    //!< HartreeOnly (pair with the Becke XC route ONLY) drops the
                                                    //!< 2/3*alpha_max XC floor so diffuse density pairs route to
                                                    //!< coarse ladder levels -- the honest diffuse-basis speedup.
};
//! The struct-parameter factory (preferred surface; the positional overload above forwards here).
Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       const GPWParams& p = {});

//! Emit the GPW `grids` report section for \a bs (the grids are k-independent, so the first block speaks for
//! all).  The orchestrator calls this AFTER the conditioning pre-flight, so a singular basis aborts before the
//! grid ladder is ever built (fail-fast; doc/RunReportPlan.md).  No-op if a run is not open.
void EmitGpwGrids(const Complex_BS& bs);

//! Fail-fast conditioning PRE-FLIGHT (doc/RunReportPlan.md, doc/GPWPlan1.md §4a): for each Bloch block of \a bs,
//! analyse its ANALYTIC, grid-free overlap and emit the report's basis.perIrrep (conditioning) + basis.removed
//! (redundant AOs).  Returns the total redundant-function count (0 == well-conditioned).  Run it under an open
//! Section("basis") BEFORE building the Hamiltonian/grids; a positive return is the cue to abort before any
//! grid work -- instead of building the whole ladder and only then discovering the basis is singular.
size_t VetGpwConditioning(const Complex_BS& bs);

} //namespace

// NOT exported: the concrete containers are an implementation detail (callers use Factory/GPWFactory's
// abstract Complex_BS).  Named here rather than anonymous in Imp/ so they are first-class types -- the home
// for the shared density-grid PeriodicGridEvaluator, GPW_BasisSet sitting beside PW_BasisSet.  Ctors in Imp/BasisSet.C.
namespace qchem::BasisSet::Lattice_3D
{

//! A tBasisSet<dcmplx> holding the plane-wave Bloch block(s); owns the IBS list (deleted with the basis).
class PW_BasisSet : public BasisSet::BasisSetImp<dcmplx>
{
public:
    //! Build one PlaneWave_IBS per Brillouin-zone k-point of \a lat at cutoff \a Ecut (a single Gamma block
    //! for an N=(1,1,1) mesh).  The basis ctor is the sole place that enumerates k -- the framework's
    //! per-irrep loop then IS the BZ sum \f$\sum_k w_k\f$.
    PW_BasisSet(const ::qchem::Lattice_3D& lat, double Ecut);
};

//! A tBasisSet<dcmplx> holding ONE \f$\Gamma\f$ GPW block (periodic Gaussians), built from a molecular
//! Gaussian basis over the cell's atoms + a density-collocation cutoff.  The GPW sibling of PW_BasisSet.
class GPW_BasisSet : public BasisSet::BasisSetImp<dcmplx>
{
public:
    GPW_BasisSet(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                 const GPWParams& p);
    //! The crystal reciprocal point group as \f${U|\tau}\f$ ops -- exposed so the composite density star-averages
    //! ρ̃ with the glide phase (IBZ).  Empty unless reduceBZ (the fold and the density symmetrization are one
    //! package).  (This is the FULL point group with fractional translations, distinct from the linear-only ops
    //! the ctor uses to fold the mesh; doc/GPWPlan1.md items 3 + 5.)
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetReciprocalPointOps() const override {return itsReciprocalOps;}
    //! The full DETECTED \f${U|\tau}\f$ set -- populated imposed or not (detection always runs; it is the
    //! §3 order-parameter diagnostic's reference group, doc/SymmetryUpgradePlan.md).
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetDetectedReciprocalOps() const override {return itsDetectedOps;}
private:
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> itsReciprocalOps;   //!< {} unless reduceBZ (doc/GPWPlan1.md items 3+5)
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> itsDetectedOps;     //!< the full detected group, always (§3 diagnostic)
};

} //namespace
