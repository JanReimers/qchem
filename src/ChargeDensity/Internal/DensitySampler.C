// File: ChargeDensity/Internal/DensitySampler.C  The two concrete sampling strategies behind
// `qchem.ChargeDensity.DensitySampler`, plus the helper both implementation units share.
//
// INTERNAL, deliberately: a client asks `MakeDensitySampler` which strategy fits its fit basis and holds the
// abstract face (user, 2026-09-10: *"Give the clients what they need and nothing more"*).  ⚠ Unit tests DO
// name these directly and that is sanctioned -- CLAUDE.md: *"Unit tests are allowed to cheat and import
// Internal stuff."*
//
// "Singles" vs "Pair" is real domain vocabulary, not an implementation accident: rho from a table of SINGLE
// basis functions Phi_gi, versus rho collocated from the density matrix through the orbital-PAIR 3-centre
// tensor.
module;
#include <cassert>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <iosfwd>
#include <map>
#include <optional>   // GroupResidual (R2.1)
#include <memory>
#include <string>
#include <functional>   // GaussianFitSpace::Evaluator
#include <vector>   // SinglesDensitySampler sigmas/flipFixed (Shubnikov S3)
export module qchem.ChargeDensity.Internal.DensitySampler;
export import qchem.ChargeDensity.DensitySampler;   // the abstract face these realise
import qchem.BasisSet.Fit_Types;
import qchem.BasisSet.Orbital_DFT_IBS;      // the fit-basis faces + FitQuadrature
import qchem.BasisSet.Orbital_1E_IBS;       // Orbital_1E_IBS<double>: the periodic A1 Gamma block (GaussianFitSpace)
import qchem.BasisSet.G_FieldEvaluator;     // G_RasterTransform -- the pair route asks its raster for size/quadrature
import qchem.Parallel;                       // LazyLock (lazy per-block caches filled from outer-loop workers)
import qchem.Fitting.FunctionFitter;        // FunctionFitter_Scalar / ScalarProjector
import qchem.ChargeDensity.Types;           // tobs_t/cobs_t/robs_t -- this library has its OWN, identical to the
                                            // qcHamiltonian typedefs the engine used to borrow (R1.0e, 2026-09-10)
                                            // with no Hamiltonian dependency of its own; it moves with the engine
import qchem.ChargeDensity;
import qchem.Mesh;                          // qcMesh::Mesh/MeshParams (the quadrature the engine integrates on)
import qchem.Mesh.SiteMultipoles;           // MultipoleFitSpace: the site multipole transforms
import qchem.Symmetry.Lattice_3D.Fold;      // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)
import qchem.Symmetry;                      // sym_t / SymMap
import qchem.SlotMap;                       // InstanceSlotMap: the per-block adjoint views, by the block's instance
import qchem.BasisSet.GradOverlap3C;   // the δ basis's gradient 3-centre overlap (GGA stage 2)
import qchem.BasisSet.PairProjection;  // ScreenedMatrixIntegrator -- the pair route's MatrixAdjoint view
import qchem.BasisSet.GSpaceDensityProjector;
import qchem.BasisSet.PointSetPullback;   // the block face the Becke mesh's pullback asks (READY 4)   // the BALL route's ρ̃ (V-ONE-OF-N 3c)
export import qchem.Mesh.Integrator;        // qcMesh::MatrixAdjoint -- the ONE face this engine names
import qchem.Types;

export namespace qchem::ChargeDensity
{

// The density names this engine consumes, pulled in EXPLICITLY rather than inherited from
// qchem.Hamiltonian.  The engine must not import the term face -- that is the dependency the extraction
// exists to break -- so it states for itself what it takes.  This list IS the engine's coupling to
// qcChargeDensity, and it is what makes qcFitting an illegal home (see the header).

// NB `ReportGridCharge()` and the `using ChargeDensity::...` pulls live in the PUBLIC module -- this unit
// re-exports it, so they are in scope here without being declared twice.

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
//! by \c MakeDensitySampler from the fit basis's capabilities, once, and LATCHED for the run (switching
//! mid-SCF would change the truncated operator, i.e. the functional).
//! \brief THE SINGLES STRATEGY: \f$\rho\f$ and \f$H_{xc}\f$ both contracted through a cached basis table
//! \f$\Phi_{gi}=\chi_i(r_g)\f$ -- the implementation that works on ANY point set, and therefore the only
//! one an atom-centred (Becke) mesh can use (doc/GPWPlan1.md "Becke XC grid").
//!
//! \f$\rho(r)\f$ comes from the density's own \f$D\f$ GEMMed against the table (no FFT, no fit --
//! pointwise \f$\rho_{DM}\ge0\f$ for aufbau D, so the \f$\rho>0\f$ guard is inert), and
//! \f$\langle i|v_{xc}|j\rangle = \Phi^\dagger\,\mathrm{diag}(w\,v_{xc})\,\Phi\f$ is its exact adjoint --
//! the two faces of \c DensitySampler over ONE table, which is what makes a mismatch unrepresentable here.
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
//! Why the δ quadrature gives no force (§13.2): the ONE wording its refusals share (MissingGradients() prints it).
inline std::string NoPullbackWhy(const std::string& block)
{return "the delta XC quadrature has no geometry gradient on the block "+block+" (it needs the block's jets at its points, and an "
        "atom-centred mesh's points and weights move with the atoms)";}
//! \note TEMPLATED on the density / block scalar (V-XC-ONE-TERM, 2026-10-31): a molecule's real δ basis samples the
//! same way a crystal's Bloch one does.  The periodic \c SinglesDensitySampler below adds the response and real-TRIM faces.
template <class T> class tSinglesDensitySampler
    : public virtual tDensitySampler<T>
{
public:
    //! \brief Built ON the \f$\delta\f$ fit basis, which IS the quadrature: it owns the points, the
    //! weights, the orbit fold and the Shubnikov tags, and ANSWERS with per-FUNCTION integrals,
    //! projections and symmetrizations rather than handing any of that out.  What is left here is pure
    //! POLICY -- which
    //! source \f$\rho\f$ comes from this iteration, per-serial caching, the spin channels, the
    //! DM-source damping -- none of which is basis business.
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS<T>> fit_t;
    //! \a quad is the SAME \c BasisSet::FitQuadrature the \f$\delta\f$ fit basis was built over, injected
    //! by the factory that created both (never taken from the basis -- it has no getter).  Two of its
    //! fields are consumed here and neither is a fitting question:
    //!  - \c mesh carries the ATOMIC site blocks -- a general-purpose observable, used here because this is
    //!    where the field it partitions is already sampled (\c SiteIntegrals);
    //!  - \c fold + \c sigmas + \c flipFixed are the crystal orbit partition and Shubnikov spin tags, which
    //!    star-average \f$\rho\f$ (and the \f$(\rho,m)\f$ pair) on every iteration.  Those reached this
    //!    strategy as \c FIT_SF_ABS::Symmetrize / \c SymmetrizeSpin until 2026-08-24 -- two members on a fit
    //!    face whose only contribution was the geometry the basis happened to own (user).  Injecting the
    //!    sibling fields the same way as the mesh removed both.
    //! Empty (default-constructed) => a free run with no partition: no star-average, and \c SiteIntegrals
    //! answers empty -- exactly as a raster quadrature does.
    tSinglesDensitySampler(fit_t, BasisSet::FitQuadrature quad={});
    double Integrate(const rvec_t& f) const override;
    size_t NumPoints() const override;
    //! My mesh's sampling slot (the δ basis samples a field ON that mesh: \c DeltaFit_IBS::Overlap -> \c SampleOn).
    SamplingTarget* FourierTarget() const override {return itsQuad.GetMesh() ? itsQuad.GetMesh()->FourierTarget() : nullptr;}
    //! \f$\rho(r_g)\f$ for \a cd's current serial (cached across the pair; rebuilt on a new serial),
    //! STAR-AVERAGED over the fold's orbits when one was supplied (exact projector on the invariant
    //! mesh -- §6a W1.  The E/H pair needs nothing else: on orbit-symmetric weights the projector is
    //! self-adjoint and \f$v(\rho_\mathrm{sym})\f$ is already symmetric, so \c Matrix below is the
    //! exact derivative untouched).
    //! ONE overload, since 2026-08-22: the density asks ME for each of its own blocks' tables (typed per
    //! block -- 3c-3), so there is no first-pass gap and no "ensure this block first" hint to pass.
    const rvec_t& Rho(const tChargeDensity<T>* cd) const override;
    //! \brief Spin channel \f$\rho_\sigma(r_g)\f$ for \a cd's current serial -- the SPIN-NATIVE sibling of
    //! \c Rho (SymmetryUpgradePlan §4 tier 4b), cached as the {↑,↓} PAIR under ONE serial (a polarized
    //! density's \c Version() forwards to its Up child, so a single scalar cache would alias the channels).
    //! A polarized composite answers per channel (its views); a spin-agnostic density (the seed) collapses to
    //! \f$\rho_\uparrow=\rho_\downarrow=\rho/2\f$ (the HalfDensity rule -- \f$v^\sigma(\tfrac\rho2,\tfrac\rho2)
    //! =v^P(\rho)\f$).  Fold star-average applies per channel (collinear: the spatial ops act channel-wise).
    const rvec_t& RhoPol(const tChargeDensity<T>* cd, const Spin& s) const override;
    XCRhoSource RhoSource() const override {return itsSource;}
    qcMesh::ImposedSymmetry Imposes() const override {return itsFit->Imposes();}   //!< the \f$\delta\f$ basis's ({e})
    //! R2.1: measured at the two projector calls (\c Rho, \c RhoPol), on the density just sampled -- the cost is one copy of the
    //! raster per refresh.  0 on a free run.
    std::optional<double> GroupResidual() const override {return itsGroupResidual;}
    //! \f$\langle i|v|j\rangle=\sum_g \overline{\Phi_{gi}}\,w_g v_g\,\Phi_{gj}\f$ over the cached table.
    hmat_t<T> Matrix(const tobs_t<T>* bs, const rvec_t& v) const override;
    //! \brief The derivatives beside \c Matrix on a FREE atom-centred mesh (lane B 2026-11-19, READY 4; the controller-approved split): the
    //! block's functions moving and my points riding their owners (\c qcMesh::Mesh::OwnerOffsets), my WEIGHTS HELD -- their motion is
    //! \c MakeQuadraturePullback's.  Per block, \f$c_g=w_gv_g\f$ through the block's \c BasisSet::PointSetPullback: the owner of each point gains
    //! its per-point term; the virial is the moment \f$-\sum_g\f$ perPoint\f$_g\otimes u_g\f$ (the point is \f$X_{owner}+u_g\f$, u not strained).
    //! REFUSED by name on an IMPOSED mesh (a fold: ρ is star-averaged on orbit-rotated grids), on a mesh without owners or a weight response,
    //! and on a block without the face.
    bool CanPullbackMatrix(const tobs_t<T>* bs) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const tobs_t<T>* bs, const rvec_t& v, const hmat_t<T>& D, const Structure* st) const override;
    bool CanVirialMatrix(const tobs_t<T>* bs) const override;
    Outcome<rmat3d_t,AD::NoGradient> MakeMatrixVirial(const tobs_t<T>* bs, const rvec_t& v, const hmat_t<T>& D, const Structure* st) const override;
    //! My weights' response: the mesh builder's (\c qcMesh::Mesh::WeightResponse), on a free mesh.
    bool HasQuadratureResponse() const override;
    Outcome<typename tDensitySampler<T>::QuadraturePullback,AD::NoGradient> MakeQuadraturePullback(const rvec_t& f, const Structure* st) const override;
    //! \f$\int w_A f\f$ per site over the INJECTED quadrature's mesh; empty when none was injected (or it
    //! has no site blocks -- a uniform grid has no atomic basins).  Ask, do not assume.  The MnO campaign's
    //! point probe (\f$m(r)\f$ 0.7 bohr off the nucleus, a spin DENSITY along one direction of an anisotropic
    //! d shell) is what the integrated moment built on this replaced -- doc/OpenWork.md Step 0a.
    //! PARTITION CAVEAT: Becke fuzzy basins are a CHOICE; the partition-free definition is R. F. W. Bader's
    //! QTAIM zero-flux basin, a wanted future feature.  The consumer reports which partition produced it.
    rvec_t SiteIntegrals(const rvec_t& f) const override;
protected:   // the periodic subclass contracts through the same fitter and tables
    //! Empty when my mesh can give the derivatives (free, owners, a weight response); otherwise why not.
    std::string MeshDerivativeGap() const;
    //! The block-generic bodies of the pullback and virial (U the block scalar): \a virial picks which.
    template <class U> Outcome<typename BasisSet::PointSetPullback<U>::Result,AD::NoGradient>
        PointSetOf(const tobs_t<U>* bs, const rvec_t& v, const hmat_t<U>& D, const Structure* st) const;
    template <class U> Outcome<rvec3vec_t,AD::NoGradient> PullbackT(const tobs_t<U>* bs, const rvec_t& v, const hmat_t<U>& D, const Structure* st) const;
    template <class U> Outcome<rmat3d_t,AD::NoGradient>   VirialT  (const tobs_t<U>* bs, const rvec_t& v, const hmat_t<U>& D, const Structure* st) const;
    template <class U> bool CanT(const tobs_t<U>* bs) const;
    // ⛔ Symmetrize / SymmetrizeSpin ARE GONE FROM THIS CLASS (2026-09-09, user).  They were two members
    // on the FITTING interface until 2026-08-24, moved here, then became one-line forwards when
    // qcMesh::FoldedMesh grew them as members (R1.0k).  A forwarder that adds nothing is not a
    // responsibility, it is a residue of where the code used to live -- and this class holds the
    // FoldedMesh (itsQuad), so the two call sites in Rho/RhoPol simply say itsQuad.Symmetrize(...).
    // ▶ THE RULE THAT FALLS OUT: symmetry projection belongs to whoever holds the FoldedMesh.  Any future
    // client wanting it asks its own FoldedMesh; nobody should re-export it.
    //! \brief MY FITTER'S PROJECTION FACE -- what a density projects itself onto (2026-08-24).
    //! The fitter holds the fit basis's \f$\Phi\f$ handles, so the \f$\rho\f$ FORWARD and the
    //! \f$H_{xc}\f$ ADJOINT come off one object instead of two callers asking the basis separately.
    //! On EITHER fit axis (the real one since 2026-11-12); null only for a fitter without the forward, which then samples pointwise.
    const BasisSet::PairProjector<T>* MatrixProjector() const;
    //! The FIELD projection face every δ fitter carries, on either fit axis -- the pointwise route.
    const Fitting::FieldProjector& FieldProjection() const;
    //! \f$\rho\f$ of a density-matrix density over my points: contracted through the fitter's matrix-forward handles when it has
    //! them (every δ fitter since 2026-11-12), else evaluated pointwise as a field.  Either way it is \f$\rho[D]\f$ (pin 28), equal
    //! to roundoff (\c UTChargeDensity's \c SinglesSampler gates).  CAPABILITY decides, never T.
    rvec_t ProjectDM(const tDM_CD<T>& dm) const;

    // R2.9(i): the four accessors above are CONST and everything they touch is a lazily-built cache, so the
    // caches are `mutable` -- the same idiom every other cache in this module already uses (tHT_Common::
    // itsCache, tDynamic_HT_Imp::itsCacheVersion, Dynamic_HF_HT_Imp::itsJKs).  Previously they were non-const
    // methods reached from const term methods through a non-const shared_ptr, which laundered the constness
    // without ever stating it.  itsFit is NOT mutable: it is construction-time and must not move.
    fit_t itsFit;                                 //!< the δ basis: my functions, their metric, their 3-centre overlap
    //! The SAME quadrature bundle the \f$\delta\f$ basis was built over -- injected, immutable, possibly
    //! empty.  Its mesh's point order IS the fit basis's function order (one object, handed to both), which
    //! is what makes \c SiteIntegrals' indexing and the fold's orbit indexing correct by construction; the
    //! asserts pin it anyway.
    BasisSet::FitQuadrature itsQuad;
    //! The δ SCALAR FITTER over that basis, from the same \c Fitting::Factory the molecular XC term uses
    //! (R1.0 Liskov conformance, 2026-08-22).  H_xc is now "fit the field, contract against this block",
    //! the same two calls on either representation -- so this strategy no longer performs a quadrature
    //! itself, it composes one.
    std::unique_ptr<Fitting::DeltaScalarFit<T>> itsScalarFitter;   //!< the PROJECTION face only (holds no field)
    //! The FIELD fits \c MatrixT contracts against: one fitter per distinct v (an iteration asks with v_up and v_dn), so
    //! concurrent Fock blocks (OOD V1.41 / doc/KAxisPlan.md P2) never refit a fitter another block is reading.  A slot is
    //! recycled only when no block holds its fitter; v is compared exactly (the fit is re-done only when v CHANGES).
    struct FieldFit { rvec_t v; std::shared_ptr<Fitting::DeltaScalarFit<T>> fitter; };
    mutable std::vector<FieldFit> itsFieldFits;
    qchem::LazyLock               itsFieldFitsLock;
    static constexpr size_t       kFieldFits=2;   //!< the two spin channels of one iteration
    std::shared_ptr<const Fitting::DeltaScalarFit<T>> FitFor(const rvec_t& v) const;
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
    //! The SEQUENCE of ρ SOURCES this engine has served XC from (\c NoteRoute), one entry per change: "DM"
    //! (a density matrix), "cusp-deficit" (N4: exact ρ[D] + the band-limited mixed difference), "DM-source"
    //! (the wholesale QCHEM_XC_DM_SOURCE route), "matrix-free" (a projected field).  Emitted into the run
    //! report as `scf.xcRhoRoutes` so a gate can PROVE which routes a recipe exercised -- a sequence, not the
    //! last value, because every run ends on "DM" (the converged density is sampled once more) and the
    //! route under test lives mid-trajectory.  The DM-source routes had zero coverage until 2026-09-15
    //! precisely because nothing could see them (doc/OpenWork.md TE, V1.18e).
    mutable std::vector<std::string> itsRoutes;
    //! Notes \a route in \c itsRoutes AND sets \c itsSource -- one call, so the two cannot disagree.
    void NoteRoute(const char* route, XCRhoSource src) const;
    mutable XCRhoSource itsSource=XCRhoSource::Fitted;   //!< the source of the cached rasters (\c RhoSource)
    mutable double itsGroupResidual=0.0;                 //!< R2.1: the cached density vs the group, before projection
    //! \brief The DM-source RUNNING MIX (QCHEM_XC_DM_MIX), in its OWN storage.
    //! \warning It must NOT live in itsRho/itsRhoUp/itsRhoDn.  Those are written by BOTH sampling routes,
    //! and the interleaving is adversarial: per iteration the Fock build blends through the DM-SOURCE
    //! branch and the energy evaluation then OVERWRITES the same buffer through the plain DM branch with
    //! rho[D_n] -- which is precisely the density the next blend mixes with, so
    //! (1-a)rho[D_n] + a rho[D_n] = rho[D_n] and the damping silently becomes the identity.  Measured
    //! 2026-08-21: alpha=0.25 and alpha=1.0 produced bit-identical runs before this was separated out.
    mutable rvec_t itsXCMix, itsXCMixUp, itsXCMixDn;
};

//! \brief The PERIODIC singles strategy: \c tSinglesDensitySampler<dcmplx> plus the faces only a Bloch run has -- the
//! linear-response transition pair and the real TRIM block (\c DensitySampler).
class SinglesDensitySampler
    : public tSinglesDensitySampler<dcmplx>
    , public virtual DensitySampler
{
public:
    using tSinglesDensitySampler<dcmplx>::tSinglesDensitySampler;
    using tSinglesDensitySampler<dcmplx>::Matrix;
    //! δρ_σ at my points, pair by pair through the δ basis's Φ tables (\c Transition_Overlap3C) -- the same
    //! tables as Rho, no cache, no orbit mean.
    cvec_t Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const override;
    //! Its adjoint on the same tables (\a q unused: the point route samples the Bloch-q function itself).
    mat_t<dcmplx> Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t& q) const override;
    //! The REAL-BLOCK sibling (Step 3c): a real TRIM block's \f$\Phi\f$ table is real, so its quadrature
    //! GEMM runs in REAL arithmetic -- the first place the Step-3 quadrature win is actually realized.
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
    using tSinglesDensitySampler<dcmplx>::MakeMatrixPullback;
    using tSinglesDensitySampler<dcmplx>::CanPullbackMatrix;
    using tSinglesDensitySampler<dcmplx>::MakeMatrixVirial;
    using tSinglesDensitySampler<dcmplx>::CanVirialMatrix;
    //! The real TRIM block's derivatives: the Bloch block's bodies on its own scalar.
    bool CanPullbackMatrix(const robs_t* bs) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
    bool CanVirialMatrix(const robs_t* bs) const override;
    Outcome<rmat3d_t,AD::NoGradient> MakeMatrixVirial(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
};

//! \brief THE GRADIENT PROJECTOR (GGA stages 2-3): the density's ONE operation (\c tChargeDensity::Project) onto the δ basis's GRADIENT
//! 3-centre overlap, so \f$\nabla\rho\f$ at the points comes off the density exactly as \f$\rho\f$ does -- the density contracts its own D, and
//! no density face changes.  Its coefficient axis is 3N long, component-major (\c BasisSet::Integrals_GradOverlap3C).  Per block it vends the
//! forward + adjoint PAIR (\c ScreenedMatrixIntegrator over the tensor), so the gradient half of H is the exact transpose of what sampled it.
//! A density with no D (the seed) projects its FIELD: the field's own analytic gradient at the points.  The \c Forward overloads are each
//! lineage's block faces (CRTP, as \c DeltaAxisFitFaces), so the periodic projector serves the Bloch AND the real TRIM blocks.
template <class Self, class T> class GradProjectorFaces;
template <class Self> class GradProjectorFaces<Self,double>
    : public virtual BasisSet::PairProjector<double>
{
public:
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override
        {return static_cast<const Self&>(*this).template Integrator<double>(orb);}
};
template <class Self> class GradProjectorFaces<Self,dcmplx>
    : public virtual BasisSet::PairProjector<dcmplx>
{
public:
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override
        {return static_cast<const Self&>(*this).template Integrator<double>(orb);}
    const qcMesh::MatrixForward<dcmplx>& Forward(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override
        {return static_cast<const Self&>(*this).template Integrator<dcmplx>(orb);}
};

template <class T> class DeltaGradientProjector
    : public GradProjectorFaces<DeltaGradientProjector<T>,T>
{
public:
    //! \a fb must have the gradient face (THROWS \c invalid_argument otherwise); \a mesh is its quadrature (one point per function).
    DeltaGradientProjector(std::shared_ptr<const BasisSet::FIT_SF_ABS<T>> fb, const qcMesh::Mesh& mesh);
    //! The block's forward + adjoint pair over its gradient tensor (\a U the block scalar), cached for the run (geometry-fixed).
    template <class U> const ScreenedMatrixIntegrator<U>& Integrator(const BasisSet::Orbital_DFT_IBS<U,T>& orb) const;
    //! \f$\nabla f(r_g)\f$, component-major: the field's own gradient (\c ScalarFunction::Gradient) at each point.
    Projections ProjectField(const ScalarFunction<double>& f) const override;
    size_t NumCoefficients() const override {return 3*itsPoints.size();}
private:
    template <class U> const InstanceSlotMap<ScreenedMatrixIntegrator<U>>& Cache() const
    {
        if constexpr (std::is_same_v<U,double>) return itsIntR;
        else                                    return itsIntC;
    }
    std::shared_ptr<const BasisSet::FIT_SF_ABS<T>>          itsFit;    //!< co-owned: the tensors' closures capture its tables
    rvec3vec_t                                              itsPoints;
    std::shared_ptr<const rvec_t>                           itsW3;     //!< the mesh weights, tiled once per component (\c Integrate only)
    InstanceSlotMap<ScreenedMatrixIntegrator<double>>       itsIntR;
    InstanceSlotMap<ScreenedMatrixIntegrator<dcmplx>>       itsIntC;   //!< used on the periodic lineage only
    qchem::LazyLock                                         itsIntLock;
};

//! \brief THE SINGLES DENSITY-GRADIENT STRATEGY (GGA stages 2-3, the point route): \a Base (the Singles strategy of the lineage) plus
//! \f$\nabla\rho_\sigma\f$ through \c DeltaGradientProjector and the \f$\{v,\mathbf{w}\}\f$ gather.  Always spin-resolved: \c Rho is the sum of the
//! channels, so one serial-keyed cache answers every face (\c tDensityGradientSampler).
//! REFUSES (throws) what it cannot sample exactly: a folded quadrature (constructor), and a density mixed as a field whose exact source is
//! retained (\c GradRhoPol: the cusp-deficit correction is a G-series whose gradient no field gives yet).
template <class T, class Base=tSinglesDensitySampler<T>> class tSinglesGradientSampler
    : public         Base
    , public virtual tDensityGradientSampler<T>
{
public:
    tSinglesGradientSampler(typename tSinglesDensitySampler<T>::fit_t fit, BasisSet::FitQuadrature quad);
    using Base::Matrix;
    //! \f$\rho=\rho_\uparrow+\rho_\downarrow\f$ off the spin cache (never the base's scalar cache: see the class note).
    const rvec_t& Rho(const tChargeDensity<T>* cd) const override;
    const rvec3vec_t& GradRhoPol(const tChargeDensity<T>* cd, const Spin& s) const override;
    hmat_t<T> Matrix(const tobs_t<T>* bs, const rvec_t& v, const rvec3vec_t& w) const override;
    //! Warms ρ AND ∇ρ, spin-resolved, whatever \a polarized says (the run's shape is the term's collapse, not this cache's).
    void WarmForDensity(const tChargeDensity<T>* cd, bool polarized) const override;
protected:
    //! The {v, w} gather on a block of scalar \a U (a Bloch block, or the periodic run's real TRIM block).
    template <class U> hmat_t<U> MatrixVW(const tobs_t<U>* bs, const rvec_t& v, const rvec3vec_t& w) const;
private:
    rvec3vec_t ProjectGradDM(const tDM_CD<T>& dm) const;
    rvec3vec_t ProjectGradField(const tChargeDensity<T>& cd) const;
    rvec3vec_t Unpack(const rvec_t& c) const;   //!< component-major 3N -> one 3-vector per point

    DeltaGradientProjector<T> itsGrad;
    mutable rvec3vec_t itsGradUp, itsGradDn;
    mutable size_t     itsGradVersion=size_t(-1);
    mutable rvec_t     itsRhoSum;
    mutable size_t     itsRhoSumVersion=size_t(-1);
};

//! \brief The PERIODIC singles density-gradient strategy (GGA stage 3): the transition and real-TRIM faces of \c SinglesDensitySampler, the
//! gradient faces above, and the real TRIM block's {v, w} gather.
class SinglesGradientSampler
    : public         tSinglesGradientSampler<dcmplx, SinglesDensitySampler>
    , public virtual DensityGradientSampler
{
public:
    using tSinglesGradientSampler<dcmplx, SinglesDensitySampler>::tSinglesGradientSampler;
    using tSinglesGradientSampler<dcmplx, SinglesDensitySampler>::Matrix;
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v, const rvec3vec_t& w) const override {return this->template MatrixVW<double>(bs, v, w);}
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
//! route and the H route were two members that happened to agree -- see the \c DensitySampler header for
//! what that permitted.
//!
//! \warning TWO ROUTES, LATCHED (R2.16).  A density-matrix-backed density projects onto the raster (\c MakeRasterDensityProjector) as
//! the RAW collocated \f$\rho_{DM}\f$; the matrix-free SEED answers EMPTY, and then \f$\rho\f$ comes from the
//! BALL round trip (\f$\tilde\rho\f$ inverse-FFT) while \f$H\f$ stays the raw adjoint.  (A plane-wave lineage, whose
//! ball holds every pair, no longer reaches this class: \c ProjectedDensitySampler, V-XC-ONE-TERM step 2, 2026-11-09.)  These minimise DIFFERENT functionals, so the route is latched
//! on the first matrix-backed density and any later change THROWS.  Iteration 0 is the one unavoidable
//! exception -- with no \f$D\f$ there is nothing to collocate -- and its energy is discarded anyway.
class PairDensitySampler
    : public virtual DensitySampler
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    //! \a fb is the raster-backed \f$v_{xc}\f$ fit basis from \c CreateVxcFitBasisSet: it supplies the
    //! quadrature (\c BasisSet::Quadrature), keys the density's collocation and the orbital's
    //! \c Overlap3C, and -- on the BALL fallback -- backs the ortho scalar fitter built here.
    explicit PairDensitySampler(fbs_t fb);
    ~PairDensitySampler();
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
    XCRhoSource RhoSource() const override {return itsSource;}
    qcMesh::ImposedSymmetry Imposes() const override {return itsFitBasis->Imposes();}   //!< the plane-wave fit on its raster (T)
    //! R2.1: the raster route does not star-average \f$\rho\f$ at all, so it cannot say how far from the group its input sat --
    //! and an imposed run on this route feeds XC an UNPROJECTED density (the Becke singles route is the one that projects).
    std::optional<double> GroupResidual() const override {return std::nullopt;}
    //! δρ_σ on the raster through the transition density's G-space face (RAW collocation, else the BALL).
    //! The periodic part \f$e^{-iq\cdot r}\delta\rho_\sigma\f$ on my raster, pair by pair through each KET block's RAW
    //! transition collocation (\c Transition_DFT_IBS::TransitionOnGrid) -- the RAW route the ground state latched;
    //! a lineage without it (a plane-wave basis) THROWS rather than falling back to the ball.
    cvec_t Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const override;
    //! Its exact transpose (\c TransitionGridAdjoint), at the wave vector \a q the field was sampled at.
    mat_t<dcmplx> Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t& q) const override;
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const override;
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
    //! The derivatives beside \c Matrix (§13.2; lane B 2026-11-19): this strategy is the RAW GPW raster gather only (the ball route is
    //! \c ProjectedDensitySampler's since V-XC-ONE-TERM step 2), so each is the block tensor's \c applyRawAdjointPullback / Virial at the
    //! fixed field -- the exact derivative of the gather \c Matrix runs (the same level fields, screener and rule).  A block whose tensor
    //! does not carry them refuses by name.  The virial is the orbital side INCLUDING the raster weight's dilation \f$\delta_{ab}\int v\rho\f$,
    //! which with the term's \f$\delta_{ab}(E-\int v\rho)\f$ (\c HasDilationVirial: the raster is fixed in fractional coordinates, weights
    //! \f$\Omega/N\f$) gives the whole fixed-D strain derivative \f$\delta_{ab}E\f$ + the gather's moment (a raw \f$\rho_g\f$ does not scale as \f$1/\Omega\f$:
    //! the Gaussians keep their shapes; it moves only with the strained centres, which the moment carries).
    bool CanPullbackMatrix(const cobs_t* bs) const override;
    bool CanPullbackMatrix(const robs_t* bs) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const cobs_t* bs, const rvec_t& v, const chmat_t& D, const Structure* st) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
    bool HasDilationVirial() const override {return true;}
    bool CanVirialMatrix(const cobs_t* bs) const override;
    bool CanVirialMatrix(const robs_t* bs) const override;
    Outcome<rmat3d_t,AD::NoGradient> MakeMatrixVirial(const cobs_t* bs, const rvec_t& v, const chmat_t& D, const Structure* st) const override;
    Outcome<rmat3d_t,AD::NoGradient> MakeMatrixVirial(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
private:
    template <class U> hmat_t<U> MatrixT(const tobs_t<U>* bs, const rvec_t& v) const;
    //! The block's raw tensor when it carries the derivative realizations (null: it refuses).
    template <class U> const PairProjection<dcmplx>* RawDerivativeTensor(const tobs_t<U>* bs) const;
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
    //! The seed's ρ̃ on this basis (iteration 0 only: a matrix-free seed has no D to collocate).  The ball ADJOINT is gone from this
    //! strategy (V-XC-ONE-TERM step 2): a plane-wave lineage runs on \c ProjectedDensitySampler.
    std::unique_ptr<BasisSet::PairProjector<dcmplx>> itsRhoProjector;
    //! The RAW ρ_DM on the same raster (R1.0n/o): \f$\langle\delta_g|\rho\rangle\f$, asked as every density answers anything,
    //! \c Project.  This strategy serves only the GPW lineage (a plane-wave one runs on \c ProjectedDensitySampler), so every
    //! D-backed block has the raw route; a matrix-free density answers its own raster values, empty for the seed.
    std::unique_ptr<BasisSet::PairProjector<dcmplx>> itsRasterProjector;
    mutable rvec_t itsRho;                        //!< ρ(r) on the raster for the current density serial
    mutable size_t itsRhoVersion=size_t(-1);      //!< density logical-clock serial itsRho was built for
    mutable bool   itsRhoIsRaw=false;             //!< is itsRho the RAW collocated ρ_DM (vs the ball round trip)?
    mutable XCRhoSource itsSource=XCRhoSource::Fitted;   //!< the source of the cached rasters (\c RhoSource)
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
    //! \brief MY ADJOINT HALF for one orbital block -- null when this lineage has no raw pair (then the
    //! BALL fallback answers).  The \c PairProjection behind it is the basis's own cached tensor; the view is
    //! cached HERE, keyed by block, because the face hands back a reference and the view borrows.
    //!
    //! ★ I HOLD THE ADJOINT AND NOTHING ELSE (2026-09-09, user: *"nobody should need both sides"*).  The
    //! FORWARD of this same tensor is what the DENSITY contracts its \f$D\f$ into; this class only needs
    //! \f$v\to\langle i|v|j\rangle\f$, so that is the only face it names.  The object is constructed with
    //! the RASTER's own energy rule, so its \c Integrate is the pinned one -- there is no second, differently
    //! ordered quadrature reachable through it.
    //! \note A template MEMBER defined in the implementation unit, like \c MatrixT beside it: both
    //! instantiations are used in that same TU, so no explicit instantiation is needed.
    template <class U> const qcMesh::MatrixAdjoint<dcmplx>*
    BlockAdjoint(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const;
    //! Per block, by the block's INSTANCE (lane H 2026-11-15: by its spatial symmetry, two blocks of equal symmetry got the first one's
    //! view).  Empty for a block with no raw pair (a plane-wave lineage), so that answer is asked once too.
    mutable InstanceSlotMap<std::optional<ScreenedMatrixIntegrator<dcmplx>>> itsAdj;
    mutable std::shared_ptr<const rvec_t> itsAdjAxis;   //!< the views' coefficient axis (zeros: the raster's rule integrates), ONE for all blocks
    qchem::LazyLock itsAdjLock;                         //!< itsAdj and itsAdjAxis are filled lazily, possibly by k-block workers
};


//! \brief THE PROJECTED STRATEGY: XC sees \f$P\rho\f$, the projection of \f$\rho[D]\f$ onto the fit space, and \f$H\f$ is its exact
//! adjoint \f$\langle ij|P\,v(P\rho)\rangle\f$ (V-XC-ONE-TERM step 2, lane S 2026-11-09; design note `lane-s-v-xc-one-term-design.md` §2).
//!
//! The \f$\rho\f$ axis of the XC product space has two values (user ruling 2026-11-09, `XCρ {Exact, Fit}`): \f$\rho[D]\f$ pointwise
//! (the singles and pair strategies) or \f$P\rho\f$ (this one).  Both are variational: \f$E=\int e(P\rho)\f$ and
//! \f$H=\partial E/\partial D=\langle ij|P\,v\rangle\f$, because \f$P\f$ is self-adjoint in the quadrature it is built on (gate:
//! \c PW_Terms.TheXCMatrixIsItsEnergysDerivative).  Where \f$P\f$ is EXACT, \f$P\rho=\rho[D]\f$ and the cell is pin-28-legal by
//! construction; elsewhere it is an explicit, non-default choice.
//!
//! REALIZED TODAY FOR ONE FIT BASIS: the plane-wave ball (\c G_RasterTransform + \c G_OrbitalPairBall).  For a plane-wave orbital basis
//! at >= 4 Ecut the ball holds every pair difference, so \f$P\f$ is exact (Parseval), and this is the route a PW run takes.  It used to
//! be the pair sampler's "BALL fallback", reached by a per-block capability test in its adjoint; it is now its own strategy, chosen
//! once by the factory, so neither strategy carries the other's route.  The A1 Gaussian fit (molecule / atom) is the next realization
//! of the same two faces (design note §4 step 3).
//!
//! \note \c Rho is the fit's coefficients of \f$\rho\f$ (the density's own projection, star-averaged) evaluated on the raster;
//! \c Matrix fits \f$v\f$ on that raster onto the ball and contracts the coefficients with the orbital basis's 3-centre overlap.
class ProjectedDensitySampler
    : public virtual DensitySampler
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    //! \a fb must be raster-backed (\c G_RasterTransform): the raster is where \f$P\rho\f$ is evaluated and integrated.
    explicit ProjectedDensitySampler(fbs_t fb);
    ~ProjectedDensitySampler();
    double Integrate(const rvec_t& f) const override;
    size_t NumPoints() const override;
    const rvec_t& Rho(const cChargeDensity* cd) const override;
    const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s) const override;
    //! \c ExactFit when \f$P\f$ is exact (the ball holds every orbital pair), \c Fitted otherwise (pin 28).
    XCRhoSource RhoSource() const override {return itsSource;}
    qcMesh::ImposedSymmetry Imposes() const override {return itsFitBasis->Imposes();}   //!< the ball (T)
    //! This route does not star-average on a quadrature (the projection does it in G space), so it cannot say how far from the group
    //! its input sat.
    std::optional<double> GroupResidual() const override {return std::nullopt;}
    //! ⛔ Not realized: a linear response on the projected route would need the transition density's own projection onto the fit
    //! (the ball), which no transition density offers today.  THROWS, naming the gap.
    cvec_t Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const override;
    //! ⛔ Its adjoint: the same gap.
    mat_t<dcmplx> Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t& q) const override;
    //! \f$\langle i|P\,v|j\rangle\f$: \f$v\f$ fitted onto the ball on my raster, contracted with the block's 3-centre overlap.
    //! A pure function of \a v: no dependency on which density was sampled last.
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const override;
    //! ⛔ A real TRIM block: the ortho scalar fitter implements FitContraction<dcmplx,dcmplx> and not <double,dcmplx>, so it THROWS (a
    //! hole in the FITTING layer, V1.34; a plane-wave run has no real TRIM blocks today).
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
    //! The derivatives beside \c Matrix (§13.2): the scalar fitter's own derivative of its contraction at the SAME fit of \a v that
    //! \c Matrix contracts (one \c DoFit, same sampling, same lock) -- a plane-wave block: an exact zero.  A real block has no
    //! contraction (V1.34), so it refuses by name.  (Lane B's ball-route override, moved here from the pair sampler, 2026-11-10.)
    bool CanPullbackMatrix(const cobs_t* bs) const override;
    bool CanPullbackMatrix(const robs_t*) const override {return false;}
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const cobs_t* bs, const rvec_t& v, const chmat_t& D, const Structure* st) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
    //! The strain derivative (CellHalf_Design §5): my raster is the fit basis's uniform FFT raster, fixed in fractional coordinates with
    //! weights Ω/N, so the quadrature side is the pure dilation; the orbital side is the scalar fitter's, at the same fit of \a v that
    //! \c Matrix contracts (a plane-wave block: an exact zero).  A real block refuses, as its pullback does.
    bool HasDilationVirial() const override {return true;}
    bool CanVirialMatrix(const cobs_t* bs) const override;
    bool CanVirialMatrix(const robs_t*) const override {return false;}
    Outcome<rmat3d_t,AD::NoGradient> MakeMatrixVirial(const cobs_t* bs, const rvec_t& v, const chmat_t& D, const Structure* st) const override;
    Outcome<rmat3d_t,AD::NoGradient> MakeMatrixVirial(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
private:
    const BasisSet::G_RasterTransform& Raster() const;
    //! \f$P\rho\f$ of ONE density object on my raster: its own projection onto the ball, then the inverse FFT.
    rvec_t SampleOne(const cChargeDensity* cd) const;
    void Refresh   (const cChargeDensity* cd) const;
    void RefreshPol(const cChargeDensity* cd) const;

    fbs_t itsFitBasis;                                          //!< the ball: fit space, raster, Overlap3C key
    bool  itsExact;                                             //!< does P hold every orbital pair (P rho = rho[D])?
    std::unique_ptr<Fitting::ScalarFit<dcmplx>> itsScalarFitter; //!< P v: the ortho (G-space) fit of a sampled field
    qchem::LazyLock itsScalarFitterLock;                         //!< Matrix refits per call: concurrent blocks take turns
    std::unique_ptr<BasisSet::PairProjector<dcmplx>> itsRhoProjector;   //!< P rho: what a density projects itself onto
    mutable rvec_t itsRho;                                      //!< P rho on the raster, for the current density serial
    mutable size_t itsRhoVersion=size_t(-1);
    //! \warning As on the other two strategies: the scalar and the spin caches have no cross-invalidation, so ONE engine answers one
    //! shape (the asserts enforce it).
    mutable rvec_t itsRhoUp, itsRhoDn;
    mutable size_t itsPolVersion=size_t(-1);
    mutable XCRhoSource itsSource=XCRhoSource::Fitted;
};


//! \brief A FIT SPACE BETWEEN AN INNER STRATEGY'S POINTS AND THE XC POINTS (V-XC-ONE-TERM, lane S 2026-11-15): the one exact linear pair that
//! carries the inner strategy's density (on its raster, or on the XC points themselves) into a fit space and out at the XC points, and its
//! adjoint.  The strategy \c FitAtPointsDensitySampler is built on.
//!   forward  \c AtPoints: rho at the inner points -> the fit's coefficients -> rho_P at the XC points;
//!   adjoint  \c OnRaster: v at the XC points (with their weights) -> the fit's coefficients -> the field at the inner points that the inner
//!            strategy's \c Matrix contracts.
//! Exactness of the pair (OnRaster = AtPoints^T with the raster weight Omega/N and the point weights folded in) is what makes H = dE/dD.
class FitSpace
{
public:
    virtual ~FitSpace() = default;
    virtual rvec_t AtPoints(const rvec_t& onRaster) const=0;
    virtual rvec_t OnRaster(const rvec_t& atPoints) const=0;
    virtual std::string Name() const=0;   //!< the report's and console's name for the fit space
    //! The fit space's own health facts, written one level in under the decorator's health scope (default: none).
    virtual void ReportHealthFacts() const {}
    //! \brief THE PIN-28 CERTIFICATE of what \c AtPoints hands XC, given the raster's own \a rasterSource: \c ExactFit only where the fit
    //! space reproduces the raster's field EXACTLY at the points (the plane-wave ball over a raster that is itself Parseval-exact in it);
    //! every other case answers \c Fitted.  Default \c Fitted: a fit space that does not answer is NOT certified (fail closed).
    virtual XCRhoSource SourceAtPoints(XCRhoSource rasterSource) const {(void)rasterSource; return XCRhoSource::Fitted;}
    //! What the fit space imposes on the potential it carries (\c DensitySamplerBase::Imposes; doc/Pins.md pin 30).  Pure: each states it.
    virtual qcMesh::ImposedSymmetry Imposes() const=0;
};

//! \brief THE PLANE-WAVE BALL as a fit space (I3, V-XC-ONE-TERM step 4): forward = the raster's forward FFT truncated to the ball and synthesized
//! at each point, rho_P(r_g) = sum_{G in ball} rho~(G) e^{iG.r_g} (\c G_FieldEvaluator::EvalField); adjoint = v~(G) = (1/Omega) sum_g w_g v_g
//! e^{-iG.r_g} (\c ProjectField), back on the raster.  rho_P is rho[D] where the ball holds every orbital pair; on a GPW block it is the
//! raster's rho[D] truncated to the ball, converging to rho[D] as the ball grows.  The cost is the direct sum, N_G x N_points each way.
class BallFitSpace : public virtual FitSpace
{
public:
    //! \a fb: a \c G_RasterTransform and a \c G_FieldEvaluator; \a points: the mesh the series is evaluated at.
    BallFitSpace(std::shared_ptr<const BasisSet::cFIT_SF_ABS> fb, BasisSet::FitQuadrature points);
    rvec_t AtPoints(const rvec_t& onRaster) const override;
    rvec_t OnRaster(const rvec_t& atPoints) const override;
    //! The ball reproduces any field inside it: rho_P = rho[D] exactly when the raster is an \c ExactFit (the ball holds every orbital pair,
    //! plane-wave orbitals).  A raster rho[D] of a Gaussian block is NOT inside any ball, so the synthesis is a truncation: \c Fitted (and
    //! \c MakeDensitySampler refuses to build it, pin 28).
    XCRhoSource SourceAtPoints(XCRhoSource rasterSource) const override
    {return rasterSource==XCRhoSource::ExactFit ? XCRhoSource::ExactFit : XCRhoSource::Fitted;}
    std::string Name() const override {return "the plane-wave ball";}
    qcMesh::ImposedSymmetry Imposes() const override {return itsFitBasis->Imposes();}   //!< the ball's (T)
private:
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> itsFitBasis;
    BasisSet::FitQuadrature itsPoints;
};

//! \brief THE PERIODIC A1 GAUSSIAN SET as a fit space (V-XC-ONE-TERM step 6, lane S 2026-11-15).  The molecular A1 Vxc set (A1_exch) on the
//! cell's sites, lattice-summed at Gamma, f_a(r) = sum_L g_a(r - R_a - L): REAL and cell-periodic, so one set serves every k (rho and v_xc
//! are Gamma-periodic).  It is a REAL Gamma block (\a a1, a GPW block over the A1 set), held here and never passed through the fit-basis
//! factory: the complex lineage's \c tBasisSet<dcmplx>::CreateVxcFitBasisSet returns a complex fit, and a real one there would change a
//! tier-1 hub (controller 2026-11-15: real values, the factory builds it, the sampler holds it).  Its metric S = <f_a|f_b>_cell is the
//! block's ANALYTIC lattice-summed overlap (controller ruling (3)); its values at the XC points come from the block's image-summed
//! evaluation, streamed in point chunks.
//!
//! THE WHOLE PAIR LIVES ON THE XC POINTS: the inner strategy is the delta basis on the same points (rho[D] exactly there), and
//!   forward   rho_P = R rho_X,   adjoint   v' = R v,   with   R x = Phi_X S^-1 Phi_X^T (w x)
//! (the inner strategy's Matrix carries the points' weights, so v' = R v exactly: R is its own adjoint in the w inner product).  First
//! built with <f_a|rho> integrated on the GPW raster instead: the raster cannot resolve A1's TIGHT functions near the cores, and S^-1 (near
//! singular) amplified the error -- F2 in a 16-bohr box (the pin-10 twin, 2026-11-15) had min rho_P -1.7 e/bohr^3 at its converged density
//! and cost 17.7 s per projection.  On the atom-centred Becke points the tight functions are resolved, and no raster is touched.
//! S^-1 drops the modes below \c kMetricRelTol x its largest eigenvalue (diffuse lattice-summed Gaussians in a small cell are nearly
//! dependent) and REPORTS a truncation (count and floor), as the radial Coulomb metric does.
class GaussianFitSpace : public virtual FitSpace
{
public:
    //! The set's values at a chunk of points (points x functions): the one thing the fit space asks of the set besides its metric.
    typedef std::function<rmat_t(const rvec3vec_t&)> Evaluator;
    //! \a a1: the real Gamma A1 block; \a points: the XC mesh, which is also the inner strategy's (the delta basis on it).
    GaussianFitSpace(std::shared_ptr<const BasisSet::Orbital_1E_IBS<double>> a1, BasisSet::FitQuadrature points);
    //! \brief The same fit space over ANY Gaussian set with an analytic metric \a S and an evaluator \a eval at points: the molecular A1
    //! fit (FittedVxc's retirement, lane S 2026-11-15) and the periodic Gamma block alike.  \a name is the report's and console's.
    //! \a imposes: what the set imposes (pin 30) -- the fit basis's own \c Imposes for the molecular / atom A1 fit.
    GaussianFitSpace(const rsmat_t& S, Evaluator eval, BasisSet::FitQuadrature points, std::string name, qcMesh::ImposedSymmetry imposes);
    //! ⛔ THROWS: a fitted rho never reaches XC (pin 28; the projected A1 cell was deleted with the ρ axis, 2026-11-15).  This fit space
    //! serves only the fitted-potential cell's v (\c OnRaster).
    rvec_t AtPoints(const rvec_t& rho) const override;
    rvec_t OnRaster(const rvec_t& v  ) const override {return Resample(v);}
    std::string Name() const override {return itsName;}
    qcMesh::ImposedSymmetry Imposes() const override {return itsImposes;}
    //! \brief R x = Phi_X S^-1 Phi_X^T (w x): \a x fitted onto the set in the XC quadrature and sampled back at the points -- the
    //! fitted-potential cell's whole operation (A1 x rho[D]: v(rho[D]) fitted for H, pin 32).
    rvec_t Resample(const rvec_t& atPoints) const;
    size_t NumFunctions() const {return itsSinv.rows();}
    //! \brief THE METRIC GAP ||S - Phi_X^T w Phi_X||_F / ||S||_F (controller 2026-11-15): R is self-adjoint in w whatever the mesh (so H = dE/dD
    //! holds), but it is a PROJECTOR (idempotent) only when the XC quadrature reproduces the analytic metric.  The gap says whether the XC mesh
    //! resolves the set's tight functions -- the defect the raster had, measured where it could recur.
    double MetricGap() const {return itsMetricGap;}
    void ReportHealthFacts() const override;
    //! The modes S^-1 drops (\c kMetricRelTol x the largest eigenvalue).  Calibration, not theory: lane B's radial Coulomb-metric floor
    //! (Radial::IrrepBasisSet, measured there on Kr-U), adopted here unmeasured for the overlap metric.
    static constexpr double kMetricRelTol=1e-12;
private:
    rvec_t PhiT (const rvec3vec_t& R, const rvec_t& f) const;   //!< Phi(R)^T f, streamed
    rvec_t Phi  (const rvec3vec_t& R, const rvec_t& c) const;   //!< Phi(R) c, streamed
    Evaluator   itsEval;
    BasisSet::FitQuadrature itsPoints;
    std::string itsName;
    qcMesh::ImposedSymmetry itsImposes;   //!< the periodic Gamma block's T, or the molecular / atom fit basis's own answer
    rsmat_t     itsSinv;
    double      itsMetricGap=0.0;
};

//! \brief THE MULTIPOLE CHANNELS as a fit space (FittedVxc's retirement, lane S 2026-11-15): the site multipole expansion on the atom-centred
//! XC mesh it was built from (\c BasisSet::MultipoleFit_IBS's transforms, without the basis): coefficients c = moments(w x) / shell weight,
//! and back at the points by the radial splines, R x = Synthesize(c).  A dual pair, not an orthogonal projection, so R is NOT self-adjoint
//! in w: it serves the fitted-potential cell (Matrix = inner.Matrix(R v), target W) and refuses the projected one (\c AtPoints throws:
//! multipole x P rho is an open X cell of VxcProductSpace).
class MultipoleFitSpace : public virtual FitSpace
{
public:
    MultipoleFitSpace(std::shared_ptr<const qcMesh::SiteMultipoles> sites, BasisSet::FitQuadrature points);
    rvec_t AtPoints(const rvec_t& rho) const override;
    rvec_t OnRaster(const rvec_t& v  ) const override {return Resample(v);}
    std::string Name() const override {return "the multipole channels";}
    qcMesh::ImposedSymmetry Imposes() const override {return qcMesh::ImposedSymmetry::None;}   //!< channels on every site: no group of their own
    rvec_t Resample(const rvec_t& atPoints) const;
private:
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    BasisSet::FitQuadrature itsPoints;
    rvec_t itsDiagonal;   //!< each channel's shell weight
};

//! \brief THE FITTED-POTENTIAL CELL on a point set (target W; V-XC-ONE-TERM step 6, and FittedVxc's retirement): XC of rho[D] at the inner
//! strategy's points, its v FITTED for H -- what \c FittedVxc did for H, as a decorator: Rho = inner.Rho (rho[D] exactly),
//! Matrix(v) = inner.Matrix(\c FitSpace::OnRaster(v)), and E = sum_g w_g e(rho[D]) by the inner quadrature (FittedVxc fitted eps for E as
//! well; that second fit is gone).  H is NOT dE/dD (doc/Pins.md pin 32): the policy layer warns on the cell and refuses what needs dE/dD
//! (\c qchem.XCRhoPolicy).  \a T: the molecule's real blocks or a crystal's Bloch blocks.
template <class T> class tFittedPotentialDensitySampler
    : public virtual tDensitySampler<T>
{
public:
    tFittedPotentialDensitySampler(std::shared_ptr<const FitSpace> fit, std::shared_ptr<const tDensitySampler<T>> inner)
        : itsFit(std::move(fit)), itsInnerT(std::move(inner))
    {
        if (!itsFit || !itsInnerT) throw std::logic_error("FittedPotentialDensitySampler: needs the fit space and its inner strategy");
    }
    double Integrate(const rvec_t& f) const override {return itsInnerT->Integrate(f);}
    size_t NumPoints() const override {return itsInnerT->NumPoints();}
    XCRhoSource RhoSource() const override {return itsInnerT->RhoSource();}
    qcMesh::ImposedSymmetry Imposes() const override {return itsFit->Imposes();}   //!< v lives in the fit space: its group
    std::optional<double> GroupResidual() const override {return itsInnerT->GroupResidual();}
    rvec_t SiteIntegrals(const rvec_t& f) const override {return itsInnerT->SiteIntegrals(f);}
    const rvec_t& Rho(const tChargeDensity<T>* cd) const override {return itsInnerT->Rho(cd);}
    const rvec_t& RhoPol(const tChargeDensity<T>* cd, const Spin& s) const override {return itsInnerT->RhoPol(cd, s);}
    hmat_t<T> Matrix(const tobs_t<T>* bs, const rvec_t& v) const override {return itsInnerT->Matrix(bs, itsFit->OnRaster(v));}
    void WarmForDensity(const tChargeDensity<T>* cd, bool polarized) const override {itsInnerT->WarmForDensity(cd, polarized);}
protected:
    std::shared_ptr<const FitSpace>             itsFit;
    std::shared_ptr<const tDensitySampler<T>>   itsInnerT;
};

//! The periodic fitted-potential cell: the real TRIM block's Matrix through the same fit, and the response faces refused by name.
class FittedPotentialDensitySampler
    : public         tFittedPotentialDensitySampler<dcmplx>
    , public virtual DensitySampler
{
public:
    FittedPotentialDensitySampler(std::shared_ptr<const FitSpace> fit, std::shared_ptr<const DensitySampler> inner);
    using tFittedPotentialDensitySampler<dcmplx>::Matrix;
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override {return itsInner->Matrix(bs, itsFit->OnRaster(v));}
    //! ⛔ Not realized on this cell.  THROW, naming the gap.
    cvec_t Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const override;
    mat_t<dcmplx> Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t& q) const override;
private:
    std::shared_ptr<const DensitySampler> itsInner;   //!< the inner strategy through its periodic face (the real block)
};

//! \brief XC OF A FITTED DENSITY AT A POINT SET (V-XC-ONE-TERM steps 4 and 6): a DECORATOR over the raster strategy the fit's raster would use
//! on its own (\a raster: the projected one where a plane-wave ball holds every orbital pair, the GPW pair collocation otherwise), so it
//! inherits that strategy's EXACT forward/adjoint pair, and adds the fit space's exact pair (\c FitSpace) on top:
//!   Rho = FitSpace::AtPoints(raster rho),  E = sum_g w_g e(rho_P(r_g)),  Matrix(v) = raster.Matrix(FitSpace::OnRaster(v)).
//! So H = dE/dD exactly (gates: \c PW_Terms.TheBallAtBeckePointsIsItsEnergysDerivative on plane-wave orbitals,
//! \c VxcProductSpace.ThePlaneWaveBallAtBeckePointsIsItsEnergysDerivative on a GPW block).  I3 learned it the hard way: first built on the
//! projected strategy's ball projector and fitter, it was exact on plane-wave orbitals and 0.7% off dE/dD on a GPW block, whose
//! projector and fitter contraction are not one pair; composing with the raster strategy's own Rho/Matrix pair is exact on both.
//! The points move with the atoms, so this route has no geometry gradient or virial yet (the defaults: it says so).
class FitAtPointsDensitySampler
    : public virtual DensitySampler
{
public:
    //! \a fit: the fit space between the raster and \a points; \a raster: the strategy over the raster; \a points: the XC mesh.
    FitAtPointsDensitySampler(std::shared_ptr<const FitSpace> fit, std::shared_ptr<const DensitySampler> raster, BasisSet::FitQuadrature points);
    double Integrate(const rvec_t& f) const override;
    size_t NumPoints() const override {return itsPoints.GetMesh()->size();}
    //! What XC sees is the fit space's synthesis, so the FIT SPACE certifies it (never the raster's own source: the raster's rho[D] truncated
    //! to a ball is not rho[D]).
    XCRhoSource RhoSource() const override {return itsFit->SourceAtPoints(itsRaster->RhoSource());}
    qcMesh::ImposedSymmetry Imposes() const override {return itsFit->Imposes();}   //!< the fit space's (the ball: T)
    std::optional<double> GroupResidual() const override {return std::nullopt;}
    const rvec_t& Rho(const cChargeDensity* cd) const override;
    const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s) const override;
    //! \f$\langle i|v|j\rangle\f$ through the exact adjoint: \a v at my points carried by the fit space onto the raster, contracted by the
    //! raster strategy.
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const override;
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
    //! The raster strategy's own warm-up (its tables), which this route samples through.
    void WarmForDensity(const cChargeDensity* cd, bool polarized) const override {itsRaster->WarmForDensity(cd, polarized);}
    //! ⛔ Not realized: a linear response on this route.  THROWS, naming the gap.
    cvec_t Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const override;
    mat_t<dcmplx> Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t& q) const override;
private:
    //! \brief THE CELL'S HEALTH, every synthesis (controller 2026-11-14): a fit space is not rho[D] -- the ball's sharp cut rings in r (pin 1),
    //! a Gaussian fit has a residual -- so the run reports, under \c xcFitAtPoints: the fit's charge error \f$\int\rho_P-N\f$ (on the points)
    //! and \f$\min\rho_P\f$ there, and WARNS when the charge error exceeds the run's target (the mesh's \c SamplingAccuracy::target, eps x
    //! atoms in its own units, applied to the charge in electrons) or \f$\min\rho_P<0\f$ (e_xc of a negative density is guarded to zero).
    void ReportHealth(const rvec_t& rhoP, double N) const;
    std::shared_ptr<const FitSpace> itsFit;
    std::shared_ptr<const DensitySampler> itsRaster;
    BasisSet::FitQuadrature itsPoints;
    //! rho_P at my points, and the raster sample each was synthesized from: the synthesis is the cost, so it reruns only when the raster
    //! strategy's sample changed (a value compare, N_raster).
    mutable rvec_t itsRho, itsRhoUp, itsRhoDn;
    mutable rvec_t itsRhoFrom, itsUpFrom, itsDnFrom;
    mutable bool   itsWarnedCharge=false, itsWarnedNegative=false, itsAnnounced=false;   //!< one console line per new condition
};

} //namespace

//======================================================================================================
// MODULE-INTERNAL HELPER -- NOT exported.  Used by BOTH implementation units of this module.
//======================================================================================================
namespace qchem::ChargeDensity
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
    virtual rvec3_t Gradient  (const rvec3_t&) const override   // a silent 0 until the hand-gradient audit (2026-11-09)
        {throw std::logic_error("SampledField is grid-bound and has no gradient");}

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

