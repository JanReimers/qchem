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
#include <vector>   // SinglesDensitySampler sigmas/flipFixed (Shubnikov S3)
export module qchem.ChargeDensity.Internal.DensitySampler;
export import qchem.ChargeDensity.DensitySampler;   // the abstract face these realise
import qchem.BasisSet.Fit_Types;
import qchem.BasisSet.Orbital_DFT_IBS;      // the fit-basis faces + FitQuadrature
import qchem.BasisSet.G_FieldEvaluator;     // G_RasterTransform -- the pair route asks its raster for size/quadrature
import qchem.Parallel;                       // LazyLock (lazy per-block caches filled from outer-loop workers)
import qchem.Fitting.FunctionFitter;        // FunctionFitter_Scalar / ScalarProjector
import qchem.ChargeDensity.Types;           // tobs_t/cobs_t/robs_t -- this library has its OWN, identical to the
                                            // qcHamiltonian typedefs the engine used to borrow (R1.0e, 2026-09-10)
                                            // with no Hamiltonian dependency of its own; it moves with the engine
import qchem.ChargeDensity;
import qchem.Mesh;                          // qcMesh::Mesh/MeshParams (the quadrature the engine integrates on)
import qchem.Symmetry.Lattice_3D.Fold;      // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)
import qchem.Symmetry;                      // sym_t / SymMap: the per-block adjoint views, by spatial symmetry
import qchem.BasisSet.PairProjection;  // ScreenedMatrixIntegrator -- the pair route's MatrixAdjoint view
import qchem.BasisSet.GSpaceDensityProjector;   // the BALL route's ρ̃ (V-ONE-OF-N 3c)
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
    //! R2.1: measured at the two projector calls (\c Rho, \c RhoPol), on the density just sampled -- the cost is one copy of the
    //! raster per refresh.  0 on a free run.
    std::optional<double> GroupResidual() const override {return itsGroupResidual;}
    //! \f$\langle i|v|j\rangle=\sum_g \overline{\Phi_{gi}}\,w_g v_g\,\Phi_{gj}\f$ over the cached table.
    hmat_t<T> Matrix(const tobs_t<T>* bs, const rvec_t& v) const override;
    //! The derivative beside \c Matrix: REFUSED by name (§13.2).  It needs the block's jets at my points, and on an atom-centred mesh
    //! the points and weights themselves move with the atoms (V-NONLINEAR-MESH).
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const tobs_t<T>* bs, const rvec_t&, const hmat_t<T>&, const Structure*) const override
    {return Outcome<rvec3vec_t,AD::NoGradient>::Fail({NoPullbackWhy(bs->Name())});}
    //! \f$\int w_A f\f$ per site over the INJECTED quadrature's mesh; empty when none was injected (or it
    //! has no site blocks -- a uniform grid has no atomic basins).  Ask, do not assume.  The MnO campaign's
    //! point probe (\f$m(r)\f$ 0.7 bohr off the nucleus, a spin DENSITY along one direction of an anisotropic
    //! d shell) is what the integrated moment built on this replaced -- doc/OpenWork.md Step 0a.
    //! PARTITION CAVEAT: Becke fuzzy basins are a CHOICE; the partition-free definition is R. F. W. Bader's
    //! QTAIM zero-flux basin, a wanted future feature.  The consumer reports which partition produced it.
    rvec_t SiteIntegrals(const rvec_t& f) const override;
protected:   // the periodic subclass contracts through the same fitter and tables
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
    const BasisSet::PairProjector<dcmplx>* MatrixProjector() const;
    //! The FIELD projection face every δ fitter carries, on either fit axis -- the pointwise route.
    const Fitting::FieldProjector& FieldProjection() const;
    //! \f$\rho\f$ of a density-matrix density over my points: contracted through the fitter's matrix-forward handles when it has
    //! them (the periodic δ fitter), else evaluated pointwise as a field -- the molecular δ route, which is exactly how
    //! \c FittedVxc sampled it.  Either way it is \f$\rho[D]\f$ (pin 28).  CAPABILITY decides, never T.
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
    //! The real block's derivative: refused, as the Bloch block's.
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const robs_t* bs, const rvec_t&, const hmat_t<double>&, const Structure*) const override
    {return Outcome<rvec3vec_t,AD::NoGradient>::Fail({NoPullbackWhy(bs->Name())});}
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
//! \warning TWO ROUTES, LATCHED (R2.16).  A density-matrix-backed density answers \c GetRhoOnGrid with
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
    //! The derivatives beside \c Matrix (§13.2): this strategy is the RAW GPW raster gather only (the ball route is
    //! \c ProjectedDensitySampler's since V-XC-ONE-TERM step 2), and its jets are not built -- so it answers \c false and refuses by
    //! name, one \c MissingGradients() entry per block (lane B's hunks, moved by lane S 2026-11-10).
    bool CanPullbackMatrix(const cobs_t*) const override {return false;}
    bool CanPullbackMatrix(const robs_t*) const override {return false;}
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const cobs_t* bs, const rvec_t& v, const chmat_t& D, const Structure* st) const override;
    Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const override;
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
    //! The seed's ρ̃ on this basis (iteration 0 only: a matrix-free seed has no D to collocate).  The ball ADJOINT is gone from this
    //! strategy (V-XC-ONE-TERM step 2): a plane-wave lineage runs on \c ProjectedDensitySampler.
    std::unique_ptr<BasisSet::PairProjector<dcmplx>> itsRhoProjector;
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
    mutable SymMap<ScreenedMatrixIntegrator<dcmplx>> itsAdj;   //!< per block, by its SPATIAL symmetry
    qchem::LazyLock itsAdjLock;                                  //!< itsAdj is filled lazily, possibly by k-block workers
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

//! \brief THE PROJECTED STRATEGY OVER AN OVERLAP-METRIC FIT (the molecular / atom A1 Gaussian basis): XC sees \f$P\rho\f$ at the XC
//! mesh's points, and \f$H\f$ is its exact adjoint (V-XC-ONE-TERM step 3, lane S 2026-11-10; the \c XCρ::Fit cell of a non-orthogonal
//! \c VxcFit).
//!
//! With \f$\Phi_{ga}=f_a(r_g)\f$ the fit functions at the mesh points, \f$w_g\f$ the mesh weights, \f$S\f$ the fit basis's analytic
//! overlap and \f$\langle ij|a\rangle\f$ its analytic 3-centre overlap:
//!  - \f$p_a=\langle f_a|\rho\rangle=\sum_{ij}D_{ij}\langle ij|a\rangle\f$ (the density projects ITSELF: \c tChargeDensity::Project),
//!    \f$c=S^{-1}p\f$, \f$P\rho(r_g)=(\Phi c)_g\f$;
//!  - \f$E=\sum_g w_g\,e(P\rho(r_g))\f$;
//!  - \f$H_{ij}=\partial E/\partial D_{ij}=\sum_a\big(S^{-1}\Phi^T(w\circ v)\big)_a\langle ij|a\rangle\f$, \f$v=e'(P\rho)\f$.
//! The pair is exact for ANY mesh: the forward and the adjoint are the same analytic tensor and the same \f$\Phi,w\f$, so the derivative
//! needs no quadrature accuracy (gate: UTHamiltonian VxcProductSpace, the A1 Pρ cell's finite difference).  What the mesh's accuracy
//! decides is how well \f$\sum_g w_g e(P\rho)\f$ approximates \f$\int e(P\rho)\f$ -- the same question as on the \f$\delta\f$ route.
//!
//! A density without a matrix (a seed, a numeric atom sum) projects its FIELD on my mesh, \f$p=\Phi^T(w\circ\rho)\f$ (\c ProjectField):
//! the iteration-0 \f$\rho\f$ only, whose energy is discarded.
//!
//! What this cell is NOT: \f$P\rho\neq\rho[D]\f$ (the A1 basis is incomplete), so it departs from pin 28 by the user's explicit choice
//! (\c xc.rho "fit").  The policy layer reports that; this class reports \c XCRhoSource::Projected.
class OverlapProjectedDensitySampler
    : public virtual tDensitySampler<double>
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;
    //! \a fb: a non-orthogonal fit basis that can evaluate its functions at points (\c Evaluatable_IBS<double>) and carries the overlap
    //! metric (\c FIT_SF_NonOrtho); THROWS naming the face it lacks.  \a quad: the XC mesh (its points and weights ARE my quadrature).
    OverlapProjectedDensitySampler(fbs_t fb, BasisSet::FitQuadrature quad);
    ~OverlapProjectedDensitySampler();
    double Integrate(const rvec_t& f) const override;
    size_t NumPoints() const override;
    const rvec_t& Rho(const rChargeDensity* cd) const override;
    const rvec_t& RhoPol(const rChargeDensity* cd, const Spin& s) const override;
    //! Always \c Projected: an A1 basis does not hold every orbital pair, so \f$P\rho\neq\rho[D]\f$, and the run asked for \f$P\rho\f$.
    XCRhoSource RhoSource() const override {return XCRhoSource::Projected;}
    std::optional<double> GroupResidual() const override {return std::nullopt;}
    //! The atom-partitioned integrals of my mesh (the XC term's site moments), as on the \f$\delta\f$ route.
    rvec_t SiteIntegrals(const rvec_t& f) const override;
    //! \f$\langle i|P\,v|j\rangle\f$: a pure function of \a v (no dependency on the density sampled last).
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
private:
    class Projector;                 //!< the \c PairProjector<double> a density projects itself onto (defined in the unit)
    rvec_t SampleOne(const rChargeDensity* cd) const;   //!< \f$\Phi S^{-1}\langle f|\rho\rangle\f$
    void Refresh   (const rChargeDensity* cd) const;
    void RefreshPol(const rChargeDensity* cd) const;

    fbs_t                                itsFitBasis;
    const BasisSet::FIT_SF_NonOrtho*     itsMetric;   //!< the same object's overlap-metric face (S^-1, the 3-centre tensor)
    BasisSet::FitQuadrature              itsQuad;     //!< the XC mesh
    rmat_t                               itsPhi;      //!< f_a(r_g): points x functions, geometry-fixed
    std::unique_ptr<Projector>           itsProjector;
    mutable rvec_t itsRho;
    mutable size_t itsRhoVersion=size_t(-1);
    //! \warning As on the other strategies: the scalar and the spin caches have no cross-invalidation, so ONE engine answers one shape.
    mutable rvec_t itsRhoUp, itsRhoDn;
    mutable size_t itsPolVersion=size_t(-1);
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

