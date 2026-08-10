// File: Hamiltonian/Internal/PWTerms.C  Plane-wave (dcmplx) Kohn-Sham Hamiltonian terms.
//
// These are the THIN terms that complete the dependency inversion: each derives from the dcmplx term
// base (cStatic_HT/cDynamic_HT in qcHamiltonian), holds the abstract orbital basis cobs_t, dynamic_casts
// it UP to the abstract BasisSet::Band_FT_IBS (G-space) capability (in qcBasisSet), and asks that high-
// level question -- "the external matrix", "the Hartree matrix for this density".  The basis owns the
// integration; the term owns no G-vectors or mesh.  Energies delegate to the density's DM_Contract.
module;
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
export module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Hamiltonian.Internal.Term;        // cStatic_HT / cDynamic_HT + their _Imp cache bases
import qchem.BasisSet.Band_FT_IBS;           // the reciprocal-space capability: Hartree/XC + external PP assembly
import qchem.BasisSet.Fit_IBS;               // cFIT_CD_ABS (the density-fit basis Vee_Hartree is built with)
import qchem.Fitting.FunctionFitter;         // FunctionFitter_Density<dcmplx> (the fitter Vee_Hartree holds, built once)
import qchem.Pseudopotential.Integrals_Pseudo;    // external-PP operator-assembly mixin + the local/separable models the term owns
import qchem.Hamiltonian.Internal.ExFunctional; // the LDA functional the XC term composes with the density
import qchem.Hamiltonian.Types;                 // cobs_t
import qchem.Structure;
import qchem.Mesh;                              // qcMesh::Mesh/MeshParams (the DeltaFittedVxc quadrature)
import qchem.Symmetry.Lattice_3D.Fold;          // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)
import qchem.Symmetry.Irrep;                    // Irrep: the Phi-table key (spatial block identity)

export namespace qchem::Hamiltonian
{

//! Process-wide diagnostic toggle (default OFF).  When true,
//! \c PWFittedVxc::RefreshRhoGrid emits a one-line report each time it (re)collocates the density: the grid-integrated
//! charge \f$\int\rho_{\text{grid}}\f$, the analytic charge \f$\mathrm{Tr}(DS)\f$, and their difference -- the
//! CHARGE LOST TO GRID TRUNCATION (== CP2K's "Electronic density on regular grids: <int> <error>" readout).
//! A cheap, controlled number for "is the density cutoff high enough" (see doc/GPWPlan.md \S0).  Flip in place:
//! `qchem::Hamiltonian::ReportGridCharge() = true;`.
bool& ReportGridCharge();

// THE LOCAL-PP RANGE SPLIT, AS THREE TERMS (naming convention: a term that carries one side of a
// short/long split SAYS SO in its name).  The local pseudopotential is split by range -- the deep-well
// erf tail (LONG) is folded through the G-space Poisson solve as a Gaussian core charge instead of a
// per-orbital-pair sharp-field sweep (the CP2K split, doc/GPWPlan.md 0e-PP), while the compact
// poly-Gaussian remainder (SHORT) rides the direct sweep.  That is a computational decomposition, so it
// used to hide inside two terms whose names claimed something else: the LONG piece lived in the Hartree
// term (making a "Hartree" term contribute to E_een) and the SHORT piece was called simply "the
// pseudopotential" (as if it were the whole PP).  Now:
//
//   Ven_PP_Short     V_loc(short)                                -> E_een   (static)
//   Ven_PP_Long      V_loc(long), the Gaussian-core-charge fold  -> E_een   (static)
//   Ven_PP_NonLocal  the KB separable projectors                 -> E_een   (static)
//   Vee_Hartree      V_H[rho_elec]                               -> E_ee    (dynamic -- the ONLY one that
//                                                                            depends on the density)
//
// The two LOCAL halves own their own dropped-G=0 alignment (E_alphaZ), Short and Long respectively; each
// is evaluated ONCE, in the ctor, so no term re-asks at run time whether its structure is periodic.  A
// configuration without a given piece simply does not ADD that term -- the term list expresses the model,
// so no term carries a runtime "do I have a pseudopotential / projectors?" test (the Ham_PP
// `if (sep) Add(PP_NonLocal)` idiom, and the molecular PP_Local/PP_NonLocal pair this now mirrors).

//! SHORT-range LOCAL (pseudo)potential term for a plane-wave basis (static, density-independent):
//! \f$V_{loc}(\text{short})\f$.  THIS is the
//! pseudo-wall: the TERM owns the pseudopotential MODEL (an abstract local form factor), and asks the
//! basis to ASSEMBLE the matrix from it (MakeLocalPotentialShort) -- physics lives Hamiltonian-side,
//! integral assembly basis-side.  The model is non-owning (the caller keeps it alive).  Pair with
//! \c Ven_PP_Long (the other half of the range split), \c Ven_PP_NonLocal (the KB projectors, when the
//! PP has them), and the kinetic/Hartree/XC terms for a full Kohn-Sham Hamiltonian.
class Ven_PP_Short
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Ven_PP_Short(const st_t& st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    st_t theStructure;
    const Pseudopotential::LocalPotential* itsLocal;   //!< local pseudopotential model (non-owning).
    double itsAlphaZ=0.0;   //!< the SHORT G=0 alignment per electron, evaluated ONCE in the ctor (0 if finite)
};

//! The KB-separable NON-LOCAL projectors of the pseudopotential (static, density-independent):
//! \f$\sum_p h_p|\beta_p\rangle\langle\beta_p|\f$.  Its own term, mirroring the molecular lineage's
//! \c PP_Local / \c PP_NonLocal pair: a purely LOCAL pseudopotential simply does not add it, rather than
//! this term carrying a "do I have projectors?" test (\c Ham_PP's `if (sep) Add(...)` idiom).
//! No G=0 alignment -- the projectors are short-ranged by construction.
class Ven_PP_NonLocal
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    //! \a nl is REQUIRED (non-owning).  A local-only pseudopotential does not construct this term at all.
    Ven_PP_NonLocal(const st_t& st, const Pseudopotential::SeparablePotential* nl);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    st_t theStructure;
    const Pseudopotential::SeparablePotential* itsSep;   //!< KB nonlocal model (non-owning).
};

//! LONG-range half of the local pseudopotential (static, density-independent): the softened-Coulomb /
//! Gaussian-core-charge matrix \f$\langle i|V_{long}|j\rangle\f$, assembled through the same
//! \c Integrals_Pseudo cross-cast \c Ven_PP_Short uses.  Electron-ion, so its energy is
//! \f$E_{een,long}=\mathrm{Tr}(D\,V_{long})\f$ with NO \f$\tfrac12\f$ (contrast the Hartree
//! double-counting factor), and it carries the LONG part's dropped-G=0 alignment.
//!
//! It is DENSITY-INDEPENDENT -- \c MakeLocalPotentialLong takes only (structure, model) -- which is why
//! it is a plain static term.  It used to be a cached side-block inside the Hartree term, summed into
//! that term's matrix and then subtracted back out of its energy; being its own term removes both the
//! fold and the subtraction.
class Ven_PP_Long
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    //! \a loc is REQUIRED (non-owning).  A run with no local PP does not construct this term at all.
    Ven_PP_Long(const st_t& st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    st_t theStructure;
    const Pseudopotential::LocalPotential* itsLocal;   //!< local pseudopotential model (non-owning).
    double itsAlphaZ=0.0;   //!< the LONG G=0 alignment per electron, evaluated ONCE in the ctor (0 if finite)
};

// The non-relativistic kinetic ENERGY term is now the T-templated Kinetic<T>
// (qchem.Hamiltonian.Internal.Kinetic); the plane-wave Hamiltonian builds Kinetic<dcmplx>.

// The ion-ion (Ewald) ENERGY term is now the T-templated IonIon<T>
// (qchem.Hamiltonian.Internal.IonIon); the plane-wave Hamiltonian builds IonIon<dcmplx>.

//! Periodic HARTREE term for a plane-wave basis (density-dependent): the classical electron-electron
//! Coulomb potential \f$V_H[\rho_{elec}]\f$ and nothing else.  The Fock block is
//! \f$\langle i|V_H|j\rangle\f$ and the energy is \f$E_{ee}=\tfrac12\mathrm{Tr}(D V_H)\f$ -- the
//! \f$\tfrac12\f$ being the electron-electron double-counting factor.  Mirrors the molecular \c FittedVee
//! in ROLE, but not in mechanism: \c FittedVee runs a charge-constrained COULOMB-METRIC (Dunlap) fit and
//! takes its energy from the robust \f$2E_{fit}-E_{fit,fit}\f$ combination, whereas here the G-space
//! projection needs no metric SOLVE (see the V1.1/V1.1b metric discussion).
//!
//! Careful with "the projection IS the fit" -- it runs together two INDEPENDENT questions:
//!   1. the METRIC: an orthonormal fit basis makes \f$S=I\f$, so the normal equations collapse to
//!      \f$c=\langle f|\rho\rangle\f$.  That is about COST and CONDITIONING, not accuracy.
//!   2. the SPAN: whether \f$\rho\f$ actually LIES in \f$\mathrm{span}\{G\}\f$.  That is what decides
//!      whether \f$\tilde\rho=\rho\f$, and orthonormality says nothing about it.
//! The answer to (2) differs by lineage.  For PLANE-WAVE orbitals \f$\rho=\psi^*\psi\f$ is exactly
//! band-limited to the difference set \f$\{G_i-G_j\}\f$, and the CD fit basis is built at the \f$4\times\f$
//! cutoff that covers it (\c PlaneWave_IBS::CreateCDFitBasisSet says so) -- so there the representation is
//! exact and no information is lost.  For GPW (GAUSSIAN orbitals) it is NOT: a Gaussian product has
//! infinite bandwidth, so a finite \f$\{G\}\f$ ball truncates it -- which is precisely what the
//! \c ReportGridCharge diagnostic measures (charge lost to grid truncation, CP2K's "Electronic density on
//! regular grids" line).  GPW's Hartree is therefore a genuine approximation, exact only as the density
//! cutoff grows.
//!
//! Neither case is a rank-2-into-rank-1 squeeze, though.  What is represented is \f$\rho(r)\f$ -- the
//! DIAGONAL \f$\rho(r,r)\f$ -- which is genuinely a function of one point.  The map \f$D\to\tilde\rho\f$ is
//! many-to-one (it sums \f$D_{ab}\f$ over each difference \f$G_b-G_a\f$), so \f$D\f$ cannot be recovered
//! from \f$\tilde\rho\f$; but Hartree and LDA only ever need the diagonal.  The full \f$\rho(r,r')\f$ WOULD
//! be needed for exact exchange -- and consistently, the periodic density NA-asserts on the HF
//! accumulators (\c IrrepCD<dcmplx>::AccumulateExchange).
//!
//! The LONG-range local-PP fold that used to live here is now its own term, \c Ven_PP_Long: it is
//! density-INDEPENDENT and electron-ION, so it belonged in neither this term's matrix nor its energy.
class Vee_Hartree
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    //! Built with the density-fit basis (from the orbital basis's factory, exactly as \c FittedVee is).
    //! No structure and no pseudopotential model: pure \f$V_H[\rho]\f$ has no use for either.
    explicit Vee_Hartree(fbs_t chargeDensityFitBasisSet);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;

    fbs_t itsFitBasis;   //!< the CD (Coulomb-metric) fit basis, handed to the density's GetRepulsion3C
};

//! Exchange-correlation term for a plane-wave basis, carrying ONE LDA functional (so a full LDA uses a
//! Dirac PWFittedVxc + a VWN PWFittedVxc, mirroring the molecular SlaterExchange+VWN split).  The matrix is the basis
//! integral of v_xc(rho(r)); the energy is integral eps_xc(rho) rho.  Both are real-space scalar fields
//! the term composes (functional o density) and hands to the basis -- the basis owns the integration.
class PWFittedVxc
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    //! Built with the Vxc fit basis obtained from the orbital basis's factory (BuildTerms creates it ONCE,
    //! never assuming orbital==fit) -- the overlap-metric sibling of Vee_Hartree.
    PWFittedVxc(const xc_t&, fbs_t vxcFitBasisSet);
    ~PWFittedVxc();
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    //! Ensure \c itsRhoGrid holds \f$\rho(r)\f$ on the fit grid for \a cd, recomputing (one inverse FFT) only
    //! on a new density serial.  Shared by MakeMatrix (fits \f$v_{xc}\f$) and GetEnergy (integrates \f$\epsilon_{xc}\rho\f$),
    //! so the transform runs ONCE per SCF iteration, whichever runs first.
    void RefreshRhoGrid(const cChargeDensity* cd) const;

    xc_t itsXc;
    fbs_t itsVxcFitBasis;   //!< the Vxc (overlap-metric) fit basis, handed to the density's GetFourierDensity
    //! The ortho scalar fitter (built once).  It OWNS the FFT quadrature grid (from the fit basis); the XC
    //! quadrature comes from the FIT basis, not the orbital basis (so relCutoff / GridCutoffFactor control it).
    //! The term borrows that ONE grid via itsScalarFitter->Grid() -- no second cross-cast of the fit basis (#7).
    std::unique_ptr<Fitting::GriddedScalarFitter> itsScalarFitter;
    mutable rvec_t itsRhoGrid;   //!< rho(r) on the fit grid for the current density (MakeMatrix builds; GetEnergy reuses)
    //! \brief Whether \c itsRhoGrid is the RAW collocated \f$\rho_{DM}\f$ (doc/GPWPlan 0.5(f2)) rather than
    //! the ball-projected round trip.  Set per refresh from the density's \c GetRhoOnGrid answer; when true,
    //! \c MakeMatrix assembles \f$H_{xc}\f$ through the raw adjoint so the E/H pair derives from the ONE raw
    //! discrete functional (and the \f$\rho>0\f$ guard becomes inert -- \f$\rho_{DM}\ge 0\f$ for aufbau D).
    mutable bool itsRhoIsRaw=false;
    //! \brief Route-stability latch (R2.16): RAW vs BALL minimise DIFFERENT functionals, so the route must
    //! not change once the SCF is running.  Latched on the first DENSITY-MATRIX-backed density (the
    //! matrix-free seed is exempt -- with no \f$D\f$ there is no \f$\rho_{DM}\f$ to collocate, so iteration
    //! 0 can only answer BALL); any later change THROWS instead of silently swapping functionals.
    mutable bool itsRouteLatched=false;
    mutable bool itsLatchedRaw=false;
};

//! Exchange-correlation term on an ATOM-CENTRED real-space quadrature (doc/GPWPlan1.md "Becke XC grid") --
//! the real-space sibling of PWFittedVxc, carrying ONE LDA functional per instance (build a Dirac + a VWN term,
//! exactly like PWFittedVxc).  The XC field is pointwise-nonlinear and sharp at the cores -- the one term the
//! uniform FFT raster serves poorly for diffuse bases (a diffuse pair x sharp field is a two-scale
//! integrand; the atom-centred mesh is dense at the cores and its point count never scales with a
//! function's diffuse reach).  Here \f$\rho(r)\f$ is evaluated ANALYTICALLY at each mesh point straight
//! off the density's ScalarFunction face (no FFT, no collocation, no fit -- pointwise \f$\rho_{DM}\ge0\f$
//! for aufbau D, so the \f$\rho>0\f$ guard is inert like the RAW route), \f$\epsilon_{xc}/v_{xc}\f$ are
//! applied per point, and \f$\langle i|v_{xc}|j\rangle = \Phi^\dagger\,\mathrm{diag}(w\,v_{xc})\,\Phi\f$
//! over the engine's cached basis table.  Hartree stays on the uniform G-space grid -- this term swaps
//! ONLY the XC quadrature.  The QUADRATURE ENGINE below is SHARED by the exchange and correlation term
//! (exactly how the PWFittedVxc pair shares one Vxc fit basis).

//! \brief The shared quadrature engine of the Becke XC pair: the mesh, the per-Bloch-block cached basis
//! tables \f$\Phi_{gi}=\chi_i(r_g)\f$ (GEOMETRY-FIXED -- built once per run per block, keyed by
//! BasisSetID), and \f$\rho\f$ at the mesh points for the current
//! density serial (built ONCE per SCF iteration for the whole pair, via cDM_CD::DM_RhoAtPoints -- the
//! density GEMMs the tables against its private \f$D\f$).  This is what makes the route O(GEMM) per
//! iteration: without it the pair re-evaluated the Bloch image sums pointwise FOUR times per iteration
//! (2 terms x (rho sample + matrix quadrature)) -- measured 4.8 s/iteration on NaF, ~all of the Becke
//! route's runtime premium.
class XC_GridEngine
{
public:
    //! The Becke route is a pure QUADRATURE -- it has no fit basis (user 2026-08-01: a zero-function
    //! pseudo-basis was a null-object smell).  The engine receives the FINISHED quadrature at
    //! construction: the mesh (already group-averaged INVARIANT on a §3-imposed run -- assembled by
    //! the Hamiltonian's Becke branch from Structure + the imposed ops) and its orbit \a fold.  An
    //! empty fold = a free run (no star-average).  Everything symmetry-shaped happens before this
    //! ctor; per iteration the engine only applies the fold's orbit-mean to ρ.
    typedef std::shared_ptr<const qcMesh::Mesh> mesh_t;
    XC_GridEngine(mesh_t, Symmetry::Lattice_3D::Fold fold = {});
    const qcMesh::Mesh& Mesh() const {return *itsMesh;}
    //! \f$\rho(r_g)\f$ for \a cd's current serial (cached across the pair; rebuilt on a new serial),
    //! STAR-AVERAGED over the fold's orbits when one was supplied (exact projector on the invariant
    //! mesh -- §6a W1.  The E/H pair needs nothing else: on orbit-symmetric weights the projector is
    //! self-adjoint and \f$v(\rho_\mathrm{sym})\f$ is already symmetric, so \c Matrix below is the
    //! exact derivative untouched).
    //! \a ensureBlock (when given) has its \f$\Phi\f$ table built FIRST, so the rho GEMM covers it even on
    //! the very first call; blocks not yet tabled self-evaluate pointwise inside the density (first pass
    //! only).  GetEnergy passes null (no basis at hand) and reuses the iteration's table.
    const rvec_t& Rho(const cChargeDensity* cd, const cobs_t* ensureBlock=nullptr) const;
    //! \brief Spin channel \f$\rho_\sigma(r_g)\f$ for \a cd's current serial -- the SPIN-NATIVE sibling of
    //! \c Rho (SymmetryUpgradePlan §4 tier 4b), cached as the {↑,↓} PAIR under ONE serial (a polarized
    //! density's \c Version() forwards to its Up child, so a single scalar cache would alias the channels).
    //! A \c cPolarized_CD answers per channel; a spin-agnostic density (the seed) collapses to
    //! \f$\rho_\uparrow=\rho_\downarrow=\rho/2\f$ (the HalfDensity rule -- \f$v^\sigma(\tfrac\rho2,\tfrac\rho2)
    //! =v^P(\rho)\f$).  Fold star-average applies per channel (collinear: the spatial ops act channel-wise).
    const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s, const cobs_t* ensureBlock=nullptr) const;
    //! \f$\langle i|v|j\rangle=\sum_g \overline{\Phi_{gi}}\,w_g v_g\,\Phi_{gj}\f$ over the cached table.
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const;
private:
    const mat_t<dcmplx>& Phi(const cobs_t* bs) const;   //!< lazily built per block (geometry-fixed)

    // R2.9(i): the four accessors above are CONST and everything they touch is a lazily-built cache, so the
    // caches are `mutable` -- the same idiom every other cache in this module already uses (tHT_Common::
    // itsCache, tDynamic_HT_Imp::itsCacheVersion, Dynamic_HF_HT_Imp::itsJKs).  Previously they were non-const
    // methods reached from const term methods through a non-const shared_ptr, which laundered the constness
    // without ever stating it.  itsMesh/itsFold are NOT mutable: they are construction-time and must not move.
    mesh_t itsMesh;                               //!< the quadrature (invariant when the fold is live)
    Symmetry::Lattice_3D::Fold itsFold;           //!< its orbit partition ({} = free run, no averaging)
    mutable std::map<Irrep,mat_t<dcmplx>> itsPhi; //!< spatial Irrep -> (npts x n) basis table
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
};

class DeltaFittedVxc
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<const XC_GridEngine> engine_t;   //!< const: all four accessors are const (R2.9(i))
    DeltaFittedVxc(const xc_t&, engine_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;
    engine_t itsEngine;   //!< the shared mesh + Phi tables + per-serial rho (one per XC pair)
};

//! SPIN-NATIVE exchange on the Becke quadrature (SymmetryUpgradePlan §4 tier 4b) -- the periodic sibling
//! of the molecular FittedVxcPol.  Exchange is CHANNEL-SEPARABLE, so one channel-native functional (a
//! spin-tagged \c SlaterExchange -- it must NOT halve \f$\rho\f$; construct with \c SlaterExchange(alpha,
//! \c Spin::Up)) serves both channels: the Fock build calls \c MakeMatrix per spin block and each fits
//! \f$v_x^\sigma=v_x(\rho_\sigma)\f$; \f$E_x=\sum_\sigma\int\epsilon_x(\rho_\sigma)\rho_\sigma\f$.
//! Shares the pair's ONE \c XC_GridEngine with the correlation term, exactly like the unpolarized pair.
class DeltaFittedVxcPol
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<ExFunctional>  xc_t;
    typedef std::shared_ptr<const XC_GridEngine> engine_t;   //!< const: all four accessors are const (R2.9(i))
    DeltaFittedVxcPol(const xc_t&, engine_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual bool          IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;       //!< channel-native (non-halving) exchange functional, shared across channels
    engine_t itsEngine;   //!< the shared mesh + Phi tables + per-serial {↑,↓} rho pair
};

//! SPIN-NATIVE correlation on the Becke quadrature -- the periodic sibling of the molecular
//! FittedVcorrPol.  Correlation does NOT separate by channel: \f$v_c^\sigma(\rho_\uparrow,\rho_\downarrow)\f$
//! couples both densities (through \f$r_s\f$ and \f$\zeta\f$), so this term evaluates the \c SpinCorrelation
//! face against BOTH channel rasters at each mesh point; \f$E_c=\int\epsilon_c(\rho_\uparrow,\rho_\downarrow)
//! (\rho_\uparrow+\rho_\downarrow)\f$.  The spin-agnostic seed collapses inside \c XC_GridEngine::RhoPol
//! (\f$\rho_\sigma=\rho/2\f$), so no term-side fallback is needed.
class DeltaFittedVcorrPol
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<SpinCorrelation> corr_t;
    typedef std::shared_ptr<const XC_GridEngine> engine_t;   //!< const: all four accessors are const (R2.9(i))
    DeltaFittedVcorrPol(const corr_t&, engine_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual bool          IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;

    corr_t   itsCorr;     //!< the spin-native correlation functional (VWN5's two-channel face)
    engine_t itsEngine;   //!< the shared mesh + Phi tables + per-serial {↑,↓} rho pair
};

} //namespace
