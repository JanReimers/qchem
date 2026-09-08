// File: Hamiltonian/Internal/XCQuadrature.C  The XC SAMPLING ENGINE: rho at a quadrature's points, and
// the exact adjoint back to a matrix.
//
// ★ THIS IS NOT HAMILTONIAN WORK, AND THE MODULE BOUNDARY NOW SAYS SO (2026-09-08).
//
// The user's definition of the Hamiltonian library (2026-09-08):
//
//     "At a very high level Hamiltonian is just: charge density in, use orbital basis and fitted functions
//      to evaluate all integrals, spit out H_ij(rho) and E(rho) for each term."
//
// By that definition a TERM is (physics) + (ask the basis for an integral) + (contract with rho).  What is
// in this file is none of those.  It is a sampling engine, and everything it owns is
// (integration GRID) x (FIT BASIS) business -- doc/Pins.md pin 2's axis pair:
//
//   * rho sampling and its per-density-serial caches (scalar, the {up,down} pair, the DM-mix buffer)
//   * the RAW-vs-BALL route latch -- a FITTING decision
//   * the Phi table contraction and the projector it goes through -- chi(r) caching, basis business
//   * orbit star-averaging of rho and the (rho,m) pair with Shubnikov spin tags -- symmetry machinery
//   * site-partitioned moments and their reporting -- an OBSERVABLE
//   * raster geometry, voxel counts, the uniform quadrature rule -- grid management
//
// ⛔ IT CANNOT MOVE TO qcFitting, AND THE EARLIER CLAIM THAT IT COULD WAS WRONG (corrected 2026-09-08).
// The engine takes a cChargeDensity in 42 places -- Version(), ProjectOnto, the cPolarized_CD and
// FourierDensity cross-casts -- and qcChargeDensity LINKS qcFitting.  Moving the engine down into
// qcFitting would therefore close a library cycle qcFitting -> qcChargeDensity -> qcFitting, which the
// linker forbids and rightly.  The earlier legality check looked only for a cycle from the Hamiltonian
// side (there is none) and for what qcFitting already links; it did not check what the ENGINE needs.
//
// ▶ THE LIBRARY HOME IS THEREFORE AN OPEN DECISION, recorded in doc/CleanupCandidates.md R1.0e.  The two
// candidates are qcChargeDensity (legal today; the rho-sampling half genuinely is its business) and a new
// leaf library between qcChargeDensity and qcHamiltonian (cleaner layering, more scaffolding).  Extracting
// the MODULE first is what makes either one a CMake change instead of surgery -- and is worth doing on its
// own account, because the boundary is now something the compiler enforces rather than a comment.
module;
#include <cassert>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>   // XC_SinglesQuadrature sigmas/flipFixed (Shubnikov S3)
export module qchem.Hamiltonian.Internal.XCQuadrature;
import qchem.BasisSet.Orbital_DFT_IBS;      // the fit-basis faces + FitQuadrature
import qchem.BasisSet.G_FieldEvaluator;     // G_RasterTransform -- the pair route asks its raster for size/quadrature
import qchem.Fitting.FunctionFitter;        // FunctionFitter_Scalar / ScalarProjector
import qchem.Hamiltonian.Types;             // cobs_t / robs_t -- a PURE TYPEDEF module over BasisSet::Orbital_1E_IBS<T>,
                                            // with no Hamiltonian dependency of its own; it moves with the engine
import qchem.ChargeDensity;
import qchem.Mesh;                          // qcMesh::Mesh/MeshParams (the quadrature the engine integrates on)
import qchem.Symmetry.Lattice_3D.Fold;      // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)
import qchem.Symmetry.Irrep;                // Irrep: the Phi-table key (spatial block identity)
import qchem.Blaze;                         // blazem::NarrowExact (the real-TRIM narrow, promoted to qcMath 2026-09-08)
import qchem.Types;

export namespace qchem::Hamiltonian
{

// The density names this engine consumes, pulled in EXPLICITLY rather than inherited from
// qchem.Hamiltonian.  The engine must not import the term face -- that is the dependency the extraction
// exists to break -- so it states for itself what it takes.  This list IS the engine's coupling to
// qcChargeDensity, and it is what makes qcFitting an illegal home (see the header).
using ChargeDensity::cChargeDensity;
using ChargeDensity::cDM_CD;

//! Process-wide diagnostic toggle (default OFF).  When true,
//! \c XC_PairQuadrature::Refresh emits a one-line report each time it (re)collocates the density: the grid-integrated
//! charge \f$\int\rho_{\text{grid}}\f$, the analytic charge \f$\mathrm{Tr}(DS)\f$, and their difference -- the
//! CHARGE LOST TO GRID TRUNCATION (== CP2K's "Electronic density on regular grids: <int> <error>" readout).
//! A cheap, controlled number for "is the density cutoff high enough" (see doc/GPWPlan.md \S0).  Flip in place:
//! `qchem::Hamiltonian::ReportGridCharge() = true;`.
bool& ReportGridCharge();

//! \brief THE XC QUADRATURE: ONE OBJECT WITH TWO ADJOINT-PAIRED FACES -- \f$\rho\f$ at the points from
//! the density, and \f$\langle i|v|j\rangle\f$ from a field at those SAME points.
//!
//! WHY THEY ARE ONE OBJECT AND NOT TWO AXES (doc/OpenWork.md, "separation of concerns in the XC terms").
//! \f$H=\partial E/\partial D\f$ holds only if the integrate-back is the EXACT ADJOINT of the collocation
//! ON THE SAME TRUNCATED OPERATOR -- what the \c GPW.RawXCConsistencyFD gate enforces.  Split the two
//! across separate owners and a mismatch becomes EXPRESSIBLE: measured 2026-08-22, an unscreened singles
//! \f$\rho\f$ paired with a screened pair \f$H\f$ sent Si from 14 to 60 iterations and moved E by 35 μHa.
//! With both faces on one type, the pairing is a class invariant instead of a convention, and "GDM just
//! works" stops resting on discipline.
//!
//! The two implementations below are a COST/TRUNCATION STRATEGY, not a semantics: on the same points both
//! compute \f$H_{ij}=\sum_g w_g v(r_g)\chi_i(r_g)\chi_j(r_g)\f$, differing only in evaluation order and
//! hence in what each truncates (the pair route's ε-screening and multigrid boxes).  Neither always wins --
//! pair scales with the SCREENED pair count, singles with \f$n_{pts}n^2\f$ -- so which one runs is decided
//! by \c MakeXCQuadrature from the fit basis's capabilities, once, and LATCHED for the run (switching
//! mid-SCF would change the truncated operator, i.e. the functional).
class XC_Quadrature
{
public:
    virtual ~XC_Quadrature() = default;
    //! \f$\int f\,d^3r\f$ for a field sampled at MY points -- the \f$E_{xc}\f$ quadrature.  A term hands
    //! back a value array and never learns where the points are (nor which kind of mesh they came from).
    //! \note POINT vocabulary is correct HERE and nowhere below it: this face IS a quadrature (its whole
    //! job is \f$\rho\f$ at points and the adjoint back), which is exactly what the fit BASIS is not.  How
    //! each strategy answers differs accordingly -- the \f$\delta\f$ one dots the coefficients with its
    //! functions' integrals (\c FIT_SF_ABS::Integrals), the raster one uses the raster's uniform rule.
    virtual double Integrate(const rvec_t& f) const=0;
    //! How many points I sample at -- for reporting only (a term's \c Write line).
    virtual size_t NumPoints() const=0;
    //! \brief \f$\rho(r_g)\f$ for \a cd's current serial, cached across the XC pair.
    //! (No "ensure this block is tabled first" hint any more: the density now asks the QUADRATURE for each
    //! of its own blocks' tables, so there is no first-pass gap for a caller to plug -- 2026-08-22.)
    virtual const rvec_t& Rho(const cChargeDensity* cd) const=0;
    //! Spin channel \f$\rho_\sigma(r_g)\f$ -- the SPIN-NATIVE sibling of \c Rho (§4 tier 4b).  Not every
    //! quadrature can answer it (the pair route has no per-spin collocation): those THROW, and the
    //! Hamiltonian's Auto rule keeps a polarized run off them.
    virtual const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s) const=0;
    //! \f$\langle i|v|j\rangle=\sum_g w_g\,\overline{\chi_i(r_g)}v_g\chi_j(r_g)\f$ -- the EXACT ADJOINT of
    //! whatever route \c Rho took, weights included (a caller passes the bare field \f$v\f$).
    virtual chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const=0;
    //! The REAL-BLOCK sibling (Step 3c): a real TRIM block's quadrature runs in REAL arithmetic.
    virtual rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const=0;
    //! \brief The per-site INTEGRATED spin moment \f$\mu_A=\int w_A(\rho_\uparrow-\rho_\downarrow)\f$, one
    //! entry per site block of my quadrature.  Default EMPTY -- a quadrature with no atomic partition (a
    //! uniform raster) has no basins to integrate over, and the caller must ask rather than assume.
    virtual rvec_t SiteMoments(const cChargeDensity*) const {return rvec_t();}
    //! \brief Pre-warm this engine's per-density caches for \a cd -- the EAGER REFRESH PHASE
    //! (doc/OpenWork.md item **KP**).  \a polarized picks which shape to warm, because the two are
    //! mutually exclusive on one engine (see the cross-invalidation warning on both implementations: an
    //! engine answers \c Rho or \c RhoPol for a run, never both).
    //!
    //! It is expressed as ONE call rather than "the term calls Rho/RhoPol itself" so that the WARMING and
    //! the SHAPE RULE stay in the engine that owns the caches; a term that reached in by name would have
    //! to know the exclusivity rule too.
    virtual void WarmForDensity(const cChargeDensity* cd, bool polarized) const
    {
        if (polarized) { RhoPol(cd, Spin::Up); RhoPol(cd, Spin::Down); }
        else           { Rho(cd); }
    }
};

//! \brief THE SINGLES STRATEGY: \f$\rho\f$ and \f$H_{xc}\f$ both contracted through a cached basis table
//! \f$\Phi_{gi}=\chi_i(r_g)\f$ -- the implementation that works on ANY point set, and therefore the only
//! one an atom-centred (Becke) mesh can use (doc/GPWPlan1.md "Becke XC grid").
//!
//! \f$\rho(r)\f$ comes from the density's own \f$D\f$ GEMMed against the table (no FFT, no fit --
//! pointwise \f$\rho_{DM}\ge0\f$ for aufbau D, so the \f$\rho>0\f$ guard is inert), and
//! \f$\langle i|v_{xc}|j\rangle = \Phi^\dagger\,\mathrm{diag}(w\,v_{xc})\,\Phi\f$ is its exact adjoint --
//! the two faces of \c XC_Quadrature over ONE table, which is what makes a mismatch unrepresentable here.
//! It wins where \f$n\f$ is small against \f$n_{pts}\f$ or where screening is weak (MnO's 4-atom cell
//! measured a Φ-sparsity ceiling of only ~2×); the pair strategy below wins where screening bites.
//!
//! \brief The shared quadrature of an XC pair: the mesh, the per-Bloch-block cached basis
//! tables \f$\Phi_{gi}=\chi_i(r_g)\f$ (GEOMETRY-FIXED -- built once per run per block, keyed by
//! BasisSetID), and \f$\rho\f$ at the mesh points for the current
//! density serial (built ONCE per SCF iteration for the whole pair, via cDM_CD::ProjectOnto -- the
//! density GEMMs the tables against its private \f$D\f$).  This is what makes the route O(GEMM) per
//! iteration: without it the pair re-evaluated the Bloch image sums pointwise FOUR times per iteration
//! (2 terms x (rho sample + matrix quadrature)) -- measured 4.8 s/iteration on NaF, ~all of the Becke
//! route's runtime premium.
class XC_SinglesQuadrature
    : public virtual XC_Quadrature
{
public:
    //! \brief Built ON the \f$\delta\f$ fit basis, which IS the quadrature: it owns the points, the
    //! weights, the orbit fold and the Shubnikov tags, and ANSWERS with per-FUNCTION integrals,
    //! projections and symmetrizations rather than handing any of that out.  What is left here is pure
    //! POLICY -- which
    //! source \f$\rho\f$ comes from this iteration, per-serial caching, the spin channels, the
    //! DM-source damping -- none of which is basis business.
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fit_t;
    //! \a quad is the SAME \c BasisSet::FitQuadrature the \f$\delta\f$ fit basis was built over, injected
    //! by the factory that created both (never taken from the basis -- it has no getter).  Two of its
    //! fields are consumed here and neither is a fitting question:
    //!  - \c mesh carries the ATOMIC site blocks -- a general-purpose observable, used here because this is
    //!    where \f$\rho_\sigma\f$ is already cached (\c SiteMoments);
    //!  - \c fold + \c sigmas + \c flipFixed are the crystal orbit partition and Shubnikov spin tags, which
    //!    star-average \f$\rho\f$ (and the \f$(\rho,m)\f$ pair) on every iteration.  Those reached this
    //!    strategy as \c FIT_SF_ABS::Symmetrize / \c SymmetrizeSpin until 2026-08-24 -- two members on a fit
    //!    face whose only contribution was the geometry the basis happened to own (user).  Injecting the
    //!    sibling fields the same way as the mesh removed both.
    //! Empty (default-constructed) => a free run with no partition: no star-average, and \c SiteMoments
    //! answers empty -- exactly as a raster quadrature does.
    XC_SinglesQuadrature(fit_t, BasisSet::FitQuadrature quad={});
    double Integrate(const rvec_t& f) const override;
    size_t NumPoints() const override;
    //! \f$\rho(r_g)\f$ for \a cd's current serial (cached across the pair; rebuilt on a new serial),
    //! STAR-AVERAGED over the fold's orbits when one was supplied (exact projector on the invariant
    //! mesh -- §6a W1.  The E/H pair needs nothing else: on orbit-symmetric weights the projector is
    //! self-adjoint and \f$v(\rho_\mathrm{sym})\f$ is already symmetric, so \c Matrix below is the
    //! exact derivative untouched).
    //! ONE overload, since 2026-08-22: the density asks ME for each of its own blocks' tables (typed per
    //! block -- 3c-3), so there is no first-pass gap and no "ensure this block first" hint to pass.
    const rvec_t& Rho(const cChargeDensity* cd) const override;
    //! \brief Spin channel \f$\rho_\sigma(r_g)\f$ for \a cd's current serial -- the SPIN-NATIVE sibling of
    //! \c Rho (SymmetryUpgradePlan §4 tier 4b), cached as the {↑,↓} PAIR under ONE serial (a polarized
    //! density's \c Version() forwards to its Up child, so a single scalar cache would alias the channels).
    //! A \c cPolarized_CD answers per channel; a spin-agnostic density (the seed) collapses to
    //! \f$\rho_\uparrow=\rho_\downarrow=\rho/2\f$ (the HalfDensity rule -- \f$v^\sigma(\tfrac\rho2,\tfrac\rho2)
    //! =v^P(\rho)\f$).  Fold star-average applies per channel (collinear: the spatial ops act channel-wise).
    const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s) const override;
    //! \f$\langle i|v|j\rangle=\sum_g \overline{\Phi_{gi}}\,w_g v_g\,\Phi_{gj}\f$ over the cached table.
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const override;
    //! The REAL-BLOCK sibling (Step 3c): a real TRIM block's \f$\Phi\f$ table is real, so its quadrature
    //! GEMM runs in REAL arithmetic -- the first place the Step-3 quadrature win is actually realized.
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
    //! \brief The per-site INTEGRATED spin moment \f$\mu_A=\int w_A(r)\,[\rho_\uparrow-\rho_\downarrow]\,d^3r\f$
    //! (electrons; \f$\times\,\mu_B\f$ for the magnetic moment), one entry per mesh site block.
    //!
    //! THE observable an atomic moment actually is — and it is FREE here: this engine already samples
    //! \f$\rho_\sigma\f$ at every mesh point once per density serial (\c RhoPol, cached), and the mesh's
    //! weights already carry the per-site partition \f$w_A\f$, so the answer is one block sum over data in
    //! hand.  It replaces the MnO campaign's point probe — \f$m(r)\f$ evaluated 0.7 bohr off the nucleus
    //! along \f$+x\f$ — which was a spin DENSITY, was never derived, and (being one direction through an
    //! anisotropic d shell) responded to the ORBITAL OCCUPATION as much as to the moment.  See
    //! doc/OpenWork.md Step 0a.
    //! \return empty when the mesh carries no site blocks (a uniform grid has no atomic partition to
    //! integrate over) — ask, do not assume.
    //! PARTITION CAVEAT: Becke fuzzy basins are a CHOICE; the partition-free definition is R. F. W. Bader's
    //! QTAIM zero-flux basin, a wanted future feature.  Report which partition produced the number.
    rvec_t SiteMoments(const cChargeDensity* cd) const override;
private:
    //! Report the current \f$\rho_\sigma\f$ pair's site moments -- called from \c RhoPol's serial-advance
    //! branch, so exactly once per NEW density and never on a cache hit.  No-op without site blocks.
    void EmitSiteMoments() const;
    //! \f$\int w_A f\f$ per site over the INJECTED quadrature's mesh; empty when none was injected (or it
    //! has no site blocks -- a uniform grid has no atomic basins).  Ask, do not assume.
    rvec_t PartitionedMoments(const rvec_t& f) const;
    //! \brief STAR-AVERAGE a coefficient vector over the injected quadrature's orbit fold, in place (§6a W1).
    //! No fold => a free run => exact no-op, so no caller asks whether symmetry was imposed.  REAL-space, so
    //! it PRESERVES \f$\rho\ge0\f$ -- XC stays on the non-negative \f$\rho_{DM}\f$ samples.
    void Symmetrize(rvec_t& f) const;
    //! \brief The MAGNETIC sibling: project the \f$(\rho,m)\f$ PAIR, which is what diagonalizes
    //! \f$\sigma\f$ -- \f$\rho\f$ EVEN under the orbit mean, \f$m\f$ ODD under the \f$\chi\f$-signed
    //! one, with the flip-fixed entries of \f$m\f$ zeroed first (Shubnikov S3, doc/SymmetryUpgradePlan.md
    //! §7).  No \f$\sigma\f$ tags => grey/free semantics => each channel averaged independently.
    void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const;
    //! \brief MY FITTER'S PROJECTION FACE -- what a density projects itself onto (2026-08-24).
    //! The fitter holds the fit basis's \f$\Phi\f$ handles, so the \f$\rho\f$ FORWARD and the
    //! \f$H_{xc}\f$ ADJOINT come off one object instead of two callers asking the basis separately.
    const Fitting::ScalarProjector& Projector() const;

    // R2.9(i): the four accessors above are CONST and everything they touch is a lazily-built cache, so the
    // caches are `mutable` -- the same idiom every other cache in this module already uses (tHT_Common::
    // itsCache, tDynamic_HT_Imp::itsCacheVersion, Dynamic_HF_HT_Imp::itsJKs).  Previously they were non-const
    // methods reached from const term methods through a non-const shared_ptr, which laundered the constness
    // without ever stating it.  itsFit is NOT mutable: it is construction-time and must not move.
    fit_t itsFit;                                 //!< the δ basis: my functions, their metric, their 3-centre overlap
    //! The SAME quadrature bundle the \f$\delta\f$ basis was built over -- injected, immutable, possibly
    //! empty.  Its mesh's point order IS the fit basis's function order (one object, handed to both), which
    //! is what makes \c SiteMoments' indexing and the fold's orbit indexing correct by construction; the
    //! asserts pin it anyway.
    BasisSet::FitQuadrature itsQuad;
    //! The δ SCALAR FITTER over that basis, from the same \c Fitting::Factory the molecular XC term uses
    //! (R1.0 Liskov conformance, 2026-08-22).  H_xc is now "fit the field, contract against this block",
    //! the same two calls on either representation -- so this strategy no longer performs a quadrature
    //! itself, it composes one.
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsScalarFitter;
    mutable rvec_t itsFittedV;                    //!< the v the fitter currently holds (refit only on change)
    //! \f$\langle f_a|1\rangle\f$ from the fit basis -- what \c Integrate dots the coefficients against.
    //! Cached because it is geometry-fixed and ~100k entries wide, and \c Integrate runs a few times per
    //! SCF iteration; REAL because a \f$\delta\f$ basis's own integrals are its (real) weights, widened to
    //! the face's Bloch scalar on the way out and narrowed once, here.
    const rvec_t& FunctionIntegrals() const;
    mutable rvec_t itsIntegrals;
    template <class U> hmat_t<U> MatrixT(const tobs_t<U>* bs, const rvec_t& v) const;
    //! \warning The scalar cache (itsRho) and the spin-resolved pair (itsRhoUp/Dn) have NO cross-
    //! invalidation: each guards only its own serial, so if one term drove \c Rho and another \c RhoPol on
    //! the SAME engine for different densities, both would report "fresh" while one held a stale raster.
    //! Unreachable today -- an engine is shared only within ONE xc/correlation PAIR and a pair is either
    //! polarized (RhoPol only) or not (Rho only) -- and the assert in each accessor pins that.  Anything
    //! that makes a run drive both (a mixed or GGA route) must add real cross-invalidation first.
    mutable rvec_t itsRho;
    mutable size_t itsRhoVersion=size_t(-1);      //!< density logical-clock serial itsRho was built for
    mutable rvec_t itsRhoUp, itsRhoDn;            //!< per-channel rasters (the polarized pair, one serial)
    mutable size_t itsPolVersion=size_t(-1);      //!< density serial the {↑,↓} pair was built for
    //! \brief STALENESS GUARD for the DM-source route.  These caches key on the MIXED FIELD's serial, but
    //! that route samples a DIFFERENT object -- the retained density matrix -- with its own serial.  If the
    //! field ever advances while its source does not, XC is served a stale rho from a cache that believes it
    //! is fresh, and the symptom (a subtly wrong V_xc) looks exactly like a degraded SCF rather than a bug.
    //! Checked LIVE rather than by \c assert, deliberately: the Step-0a site-block defect was invisible for
    //! months precisely because its guard was an assert compiled out under NDEBUG, and every benchmark row
    //! is a Release run.
    mutable size_t itsSrcVersion=size_t(-1);      //!< the DM source's serial when the pair was last built
    //! \brief The DM-source RUNNING MIX (GPW_XC_DM_MIX), in its OWN storage.
    //! \warning It must NOT live in itsRho/itsRhoUp/itsRhoDn.  Those are written by BOTH sampling routes,
    //! and the interleaving is adversarial: per iteration the Fock build blends through the DM-SOURCE
    //! branch and the energy evaluation then OVERWRITES the same buffer through the plain DM branch with
    //! rho[D_n] -- which is precisely the density the next blend mixes with, so
    //! (1-a)rho[D_n] + a rho[D_n] = rho[D_n] and the damping silently becomes the identity.  Measured
    //! 2026-08-21: alpha=0.25 and alpha=1.0 produced bit-identical runs before this was separated out.
    mutable rvec_t itsXCMix, itsXCMixUp, itsXCMixDn;
};

//! \brief THE PAIR STRATEGY: \f$\rho\f$ COLLOCATED from the density-matrix through the orbital-pair
//! 3-centre tensor, and \f$H_{xc}\f$ through that same tensor's RAW ADJOINT -- the production GPW route.
//!
//! Its points are the fit basis's FFT raster (fractional corners \f$i/N\f$, weight \f$\Omega/N_{pts}\f$),
//! taken through \c BasisSet::Quadrature exactly like any other mesh; the pair machinery needs the fit
//! BASIS as well, because \c Overlap3C(fitBasis) is what carries \c applyRaw / \c applyRawAdjoint.  That
//! is the one structural asymmetry between the two strategies: a pair quadrature takes the fit basis
//! where a singles quadrature takes only the mesh.
//!
//! It cost pair-vs-singles nothing to unify, but it EARNED the unification: the RAW forward and the RAW
//! adjoint are box-truncated per multigrid level with the same ε-screening, so \f$H_{xc}\f$ is
//! \f$\partial E_{xc}/\partial D\f$ of the ONE raw discrete functional to machine precision (gate:
//! \c GPW.RawXCConsistencyFD).  Previously the same code lived in a TERM (\c PWFittedVxc), where the ρ
//! route and the H route were two members that happened to agree -- see the \c XC_Quadrature header for
//! what that permitted.
//!
//! \warning TWO ROUTES, LATCHED (R2.16).  A density-matrix-backed density answers \c GetRhoOnGrid with
//! the RAW collocated \f$\rho_{DM}\f$; a plane-wave density or the matrix-free SEED answers EMPTY, and
//! then \f$\rho\f$ comes from the BALL round trip (\f$\tilde\rho\f$ inverse-FFT) and \f$H\f$ from the
//! ortho fitter, which is NON-variational.  These minimise DIFFERENT functionals, so the route is latched
//! on the first matrix-backed density and any later change THROWS.  Iteration 0 is the one unavoidable
//! exception -- with no \f$D\f$ there is nothing to collocate -- and its energy is discarded anyway.
class XC_PairQuadrature
    : public virtual XC_Quadrature
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    //! \a fb is the raster-backed \f$v_{xc}\f$ fit basis from \c CreateVxcFitBasisSet: it supplies the
    //! quadrature (\c BasisSet::Quadrature), keys the density's collocation and the orbital's
    //! \c Overlap3C, and -- on the BALL fallback -- backs the ortho scalar fitter built here.
    explicit XC_PairQuadrature(fbs_t fb);
    ~XC_PairQuadrature();
    double Integrate(const rvec_t& f) const override;
    size_t NumPoints() const override;
    const rvec_t& Rho(const cChargeDensity* cd) const override;
    //! \brief \f$\rho_\sigma(r)\f$ on the raster — SPIN-NATIVE since 2026-08-28.
    //!
    //! ⚠ IT USED TO THROW, and the throw was the tail of a CONFLATION (user, 2026-08-28: *"polarization
    //! and XC grids ... in my mind they have nothing to do with each other ... the user should be able to
    //! select any XC grid, and pol and unpol systems, with no if statements in the code blocking that"*).
    //! Because this route could not answer per channel, \c VxcFit::Auto sent EVERY polarized run to the
    //! δ/singles route whatever its grid — so a polarized run could never take the collocation route,
    //! which is the one CP2K uses and the only variational one.  Measured cost of that coupling on MnO
    //! with the Becke mesh vetoed: 1805 s CPU and 4.5 GB, against 584 s and 491 MB, essentially all of it
    //! Φ tables over 571787 uniform points.
    //!
    //! There was never anything spin-specific in the way: \c applyRaw takes a \f$D\f$, so a channel is
    //! one call with that channel's density, and the ADJOINT needs no change at all — \c Matrix already
    //! takes a bare field, which is exactly what \f$v_{xc,\sigma}\to H_{xc,\sigma}\f$ wants.  What was
    //! missing was the per-channel CACHE and the channel walk, both of which the singles route already had.
    const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s) const override;
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const override;
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
private:
    template <class U> hmat_t<U> MatrixT(const tobs_t<U>* bs, const rvec_t& v) const;
    //! Ensure \c itsRho holds \f$\rho(r)\f$ on the raster for \a cd, recomputing only on a new density
    //! serial -- so the XC pair's two terms and their energies share ONE collocation per iteration.
    void Refresh(const cChargeDensity* cd) const;
    //! The spin-resolved sibling: fill \c itsRhoUp / \c itsRhoDn for \a cd, once per density serial.
    void RefreshPol(const cChargeDensity* cd) const;
    //! \brief Sample ONE density object onto the raster, RAW if it can collocate and BALL otherwise.
    //! The single place that decision is made, so the scalar and the two spin channels cannot diverge.
    rvec_t SampleOne(const cChargeDensity* cd, bool& isRaw) const;
    //! Latch the RAW/BALL route on the first matrix-backed density and THROW if it ever changes (R2.16 --
    //! the two routes minimise DIFFERENT discrete functionals, so switching mid-SCF moves the target).
    void LatchRoute(const cChargeDensity* cd, bool isRaw) const;

    fbs_t itsFitBasis;      //!< the raster fit basis: quadrature, collocation key, Overlap3C key
    //! The ortho scalar fitter over that basis -- used ONLY by the BALL fallback (the RAW route fits
    //! nothing: "no ball fit anywhere").  Built once; its own grid is this basis's raster.
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsScalarFitter;
    mutable rvec_t itsRho;                        //!< ρ(r) on the raster for the current density serial
    mutable size_t itsRhoVersion=size_t(-1);      //!< density logical-clock serial itsRho was built for
    mutable bool   itsRhoIsRaw=false;             //!< is itsRho the RAW collocated ρ_DM (vs the ball round trip)?
    mutable bool   itsRouteLatched=false;         //!< has a matrix-backed density fixed the route yet?
    mutable bool   itsLatchedRaw=false;           //!< ... and to which one
    //! \warning The scalar cache (\c itsRho) and the spin pair have NO cross-invalidation -- the same
    //! warning the singles route carries, and the same asserts enforce it: ONE engine answers one shape.
    mutable rvec_t itsRhoUp, itsRhoDn;            //!< per-channel rasters (the polarized pair, one serial)
    mutable size_t itsPolVersion=size_t(-1);      //!< density serial the {↑,↓} pair was built for
    //! The fit basis's RASTER face -- where the voxel count and the uniform quadrature rule live.  Not on
    //! the fit face: a plane-wave basis counts \f$\{G\}\f$ FUNCTIONS, and its raster has more voxels than
    //! it has functions, so a caller holding a raster array must ask the raster (2026-08-23).
    const BasisSet::G_RasterTransform& Raster() const;
};

//! \brief Pick the assembly strategy for \a fb -- CAPABILITY decides, and the answer is fixed for the run.
//!
//! A δ fit basis carries points and nothing else, so it can only be contracted
//! through a Φ table: SINGLES.  A raster-backed one additionally carries the FFT transforms and keys the
//! orbital's 3-centre tensor, so the PAIR route is available and is chosen -- it is the production GPW
//! path and the one whose screening pays on large cells.  There is no extra input to supply: both routes
//! are already functions of (orbital basis, fit basis).
std::shared_ptr<const XC_Quadrature>
MakeXCQuadrature(const std::shared_ptr<const BasisSet::cFIT_SF_ABS>& fb,
                 BasisSet::FitQuadrature quad={});

} //namespace

//======================================================================================================
// MODULE-INTERNAL HELPER -- NOT exported.  Used by BOTH implementation units of this module.
//======================================================================================================
namespace qchem::Hamiltonian
{

// A field ALREADY SAMPLED at the quadrature's points, presented as the ProjectedScalar_R the ortho scalar
// fitter consumes (the BALL route's only client).  The values are v_xc(rho(r_g)), computed by the TERM --
// which is where the functional lives.
// It is GRID-BOUND: only the ortho fitter samples it, in bulk, on exactly the points it was built for.
class SampledField
    : public virtual ScalarFunction<double>
    , public         Fitting::ProjectedScalar_R
{
public:
    //! \a npts is the fit basis's own point count -- the only thing this field needs to know about the
    //! quadrature, now that the BASIS does the sampling.
    SampledField(const rvec_t& vals, size_t npts) : itsVals(vals), itsNPts(npts) {}

    // Pointwise is NOT supported: this field carries only grid values, and nothing samples it pointwise (the
    // ortho fitter uses the bulk overload).  Make the grid-bound contract explicit rather than silently wrong.
    virtual double  operator()(const rvec3_t&) const override
        {throw std::logic_error("SampledField is grid-bound: sample it in bulk on the fit grid, not pointwise");}
    virtual rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}

    // Bulk: the precomputed values, which were computed AT the fit basis's own points -- so the only thing
    // that can go wrong is a different point COUNT, and that is what the assert pins.
    virtual rvec_t  operator()(const rvec3vec_t& rs) const override
    {
        assert(rs.size()==itsNPts && itsVals.size()==itsNPts &&
               "SampledField: sampled on a different point set than the values were computed on");
        return itsVals;
    }

    virtual const ScalarFunction<double>* GetScalarFunction() const override {return this;}
private:
    const rvec_t& itsVals;   // precomputed v_xc at the fit basis's points (owned by the caller; transient)
    size_t        itsNPts;
};

} //namespace
