// File: Hamiltonian/Internal/PWTerms.C  Plane-wave (dcmplx) Kohn-Sham Hamiltonian terms.
//
// These are the THIN terms that complete the dependency inversion: each derives from the dcmplx term
// base (cStatic_HT/cDynamic_HT in qcHamiltonian), holds the abstract orbital basis cobs_t, dynamic_casts
// it UP to the abstract BasisSet::Orbital_DFT_IBS<dcmplx> (G-space) capability (in qcBasisSet), and asks that high-
// level question -- "the external matrix", "the Hartree matrix for this density".  The basis owns the
// integration; the term owns no G-vectors or mesh.  Energies delegate to the density's DM_Contract.
module;
#include <iosfwd>
#include <map>
#include <set>      // Ven_PP_NonLocal::itsByLSeen (the GPW_NL_PER_L diagnostic)
#include <vector>   // XC_SinglesQuadrature sigmas/flipFixed (Shubnikov S3)
#include <memory>
#include <string>
export module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Hamiltonian.Internal.Term;        // cStatic_HT / cDynamic_HT + their _Imp cache bases
import qchem.BasisSet.Orbital_DFT_IBS;           // the reciprocal-space capability: Hartree/XC + external PP assembly
import qchem.BasisSet.G_FieldEvaluator;      // G_RasterTransform -- the pair route asks its raster for size/quadrature
import qchem.Fitting.FunctionFitter;         // FunctionFitter_Density<dcmplx> (the fitter Vee_Hartree holds, built once)
import qchem.Pseudopotential.Integrals_Pseudo;    // external-PP operator-assembly mixin + the local/separable models the term owns
import qchem.Hamiltonian.Internal.ExFunctional; // the LDA functional the XC term composes with the density
import qchem.Hamiltonian.Types;                 // cobs_t
import qchem.Structure;
import qchem.Mesh;                              // qcMesh::Mesh/MeshParams (the XC quadrature the terms integrate on)
import qchem.Symmetry.Lattice_3D.Fold;          // Fold + SymmetrizeValues (the Becke rho star-average, §6a W1)
import qchem.Symmetry.Irrep;                    // Irrep: the Phi-table key (spatial block identity)

export namespace qchem::Hamiltonian
{

//! Process-wide diagnostic toggle (default OFF).  When true,
//! \c XC_PairQuadrature::Refresh emits a one-line report each time it (re)collocates the density: the grid-integrated
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
    , public         Static_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    Ven_PP_Short(const st_t& st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&) const;   // Step 3c: the real TRIM block
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&) const;   // the ONE assembly body
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
    , public         Static_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    //! \a nl is REQUIRED (non-owning).  A local-only pseudopotential does not construct this term at all.
    Ven_PP_NonLocal(const st_t& st, const Pseudopotential::SeparablePotential* nl);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&) const;   // Step 3c: the real TRIM block
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&) const;   // the ONE assembly body
    st_t theStructure;
    const Pseudopotential::SeparablePotential* itsSep;   //!< KB nonlocal model (non-owning).
    //! GPW_NL_PER_L=1 diagnostic (doc/SphericalLatticePlan.md I0): the per-angular-channel KB blocks,
    //! l -> (BasisSetID -> block), filled lazily by MakeMatrix and contracted per l in GetEnergy so the
    //! \f$E_{NL}^{(l)}\f$ decomposition prints beside the lumped \c EenNL.  Empty when the knob is off.
    mutable std::map<int,std::map<std::string,chmat_t>> itsByL;
    mutable std::set<std::string>                       itsByLSeen;   //!< irrep blocks already decomposed
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
    , public         Static_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    //! \a loc is REQUIRED (non-owning).  A run with no local PP does not construct this term at all.
    Ven_PP_Long(const st_t& st, const Pseudopotential::LocalPotential* loc);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&) const;   // Step 3c: the real TRIM block
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&) const;   // the ONE assembly body
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
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
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
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    fbs_t itsFitBasis;   //!< the CD (Coulomb-metric) fit basis, handed to the density's GetRepulsion3C
};

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
    //! THROWS: the pair route has no per-spin collocation (a spin-native pair quadrature is not designed).
    //! The Hamiltonian's \c VxcFit::Auto rule sends every polarized run to the δ/singles route instead.
    const rvec_t& RhoPol(const cChargeDensity*, const Spin&) const override;
    chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const override;
    rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const override;
private:
    template <class U> hmat_t<U> MatrixT(const tobs_t<U>* bs, const rvec_t& v) const;
    //! Ensure \c itsRho holds \f$\rho(r)\f$ on the raster for \a cd, recomputing only on a new density
    //! serial -- so the XC pair's two terms and their energies share ONE collocation per iteration.
    void Refresh(const cChargeDensity* cd) const;

    fbs_t itsFitBasis;      //!< the raster fit basis: quadrature, collocation key, Overlap3C key
    //! The ortho scalar fitter over that basis -- used ONLY by the BALL fallback (the RAW route fits
    //! nothing: "no ball fit anywhere").  Built once; its own grid is this basis's raster.
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsScalarFitter;
    mutable rvec_t itsRho;                        //!< ρ(r) on the raster for the current density serial
    mutable size_t itsRhoVersion=size_t(-1);      //!< density logical-clock serial itsRho was built for
    mutable bool   itsRhoIsRaw=false;             //!< is itsRho the RAW collocated ρ_DM (vs the ball round trip)?
    mutable bool   itsRouteLatched=false;         //!< has a matrix-backed density fixed the route yet?
    mutable bool   itsLatchedRaw=false;           //!< ... and to which one
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

//! \brief THE exchange-correlation term of a periodic Kohn-Sham Hamiltonian, carrying ONE LDA functional
//! (a full LDA is a Dirac instance + a VWN instance, mirroring the molecular SlaterExchange+VWN split).
//!
//! It owns the PHYSICS and nothing else: map the functional over \f$\rho\f$ at the quadrature's points,
//! hand the resulting field back for the adjoint assembly, and integrate \f$\int\epsilon_{xc}\rho\f$ on
//! the same weights.  WHICH points, WHICH representation and WHICH assembly strategy are all inside the
//! \c XC_Quadrature it was built with, so this one term serves every combination -- δ on Becke, δ on the
//! uniform cell mesh, plane-wave on the raster.
//!
//! It was \c DeltaFittedVxc, the Becke-route term, while the raster route had a term of its own
//! (\c PWFittedVxc) that duplicated this logic around its own ρ/H pair.  Two terms for one formula
//! \f$H_{ij}=\sum_g w_g v(r_g)\chi_i\chi_j\f$: the difference between them was never the physics, only the
//! evaluation order, which is exactly what \c XC_Quadrature's two implementations now hold.
class Vxc_Quadrature
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<const XC_Quadrature> quad_t;   //!< const: every accessor is const (R2.9(i))
    Vxc_Quadrature(const xc_t&, quad_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;
    quad_t   itsQuad;   //!< the shared mesh + Phi tables + per-serial rho (one per XC pair)
};

//! SPIN-NATIVE exchange on the Becke quadrature (SymmetryUpgradePlan §4 tier 4b) -- the periodic sibling
//! of the molecular FittedVxcPol.  Exchange is CHANNEL-SEPARABLE, so one channel-native functional (a
//! spin-tagged \c SlaterExchange -- it must NOT halve \f$\rho\f$; construct with \c SlaterExchange(alpha,
//! \c Spin::Up)) serves both channels: the Fock build calls \c MakeMatrix per spin block and each fits
//! \f$v_x^\sigma=v_x(\rho_\sigma)\f$; \f$E_x=\sum_\sigma\int\epsilon_x(\rho_\sigma)\rho_\sigma\f$.
//! Shares the pair's ONE \c XC_Quadrature with the correlation term, exactly like the unpolarized pair.
class Vxc_QuadraturePol
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    //! The atom-centred partition lives on my quadrature, so I am the term that can answer this
    //! (doc/OpenWork.md N1/T2).  Empty when the quadrature has no site blocks (a uniform raster).
    virtual rvec_t SiteMoments(const cChargeDensity* cd) const override;
    typedef std::shared_ptr<ExFunctional>  xc_t;
    typedef std::shared_ptr<const XC_Quadrature> quad_t;   //!< const: every accessor is const (R2.9(i))
    Vxc_QuadraturePol(const xc_t&, quad_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual bool          IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    xc_t     itsXc;       //!< channel-native (non-halving) exchange functional, shared across channels
    quad_t   itsQuad;   //!< the shared mesh + Phi tables + per-serial {↑,↓} rho pair
};

//! SPIN-NATIVE correlation on the Becke quadrature -- the periodic sibling of the molecular
//! FittedVcorrPol.  Correlation does NOT separate by channel: \f$v_c^\sigma(\rho_\uparrow,\rho_\downarrow)\f$
//! couples both densities (through \f$r_s\f$ and \f$\zeta\f$), so this term evaluates the \c SpinCorrelation
//! face against BOTH channel rasters at each mesh point; \f$E_c=\int\epsilon_c(\rho_\uparrow,\rho_\downarrow)
//! (\rho_\uparrow+\rho_\downarrow)\f$.  The spin-agnostic seed collapses inside \c XC_Quadrature::RhoPol
//! (\f$\rho_\sigma=\rho/2\f$), so no term-side fallback is needed.
class Vcorr_QuadraturePol
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp   // real TRIM block capability (Step 3c)
{
public:
    typedef std::shared_ptr<SpinCorrelation> corr_t;
    typedef std::shared_ptr<const XC_Quadrature> quad_t;   //!< const: every accessor is const (R2.9(i))
    Vcorr_QuadraturePol(const corr_t&, quad_t);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual bool          IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t MakeMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const;   // Step 3c
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    corr_t   itsCorr;     //!< the spin-native correlation functional (VWN5's two-channel face)
    quad_t   itsQuad;   //!< the shared mesh + Phi tables + per-serial {↑,↓} rho pair
};


//! \brief The exchange+correlation TERM PAIR for a run whose \f$v_{xc}\f$ fit basis is \a fb -- ready to
//! \c Add (ownership passes with each release()).
//!
//! A Hamiltonian builder asks for XC terms and gets XC terms.  WHICH quadrature they run on, which
//! assembly strategy that quadrature uses, and the fact that the pair SHARES one of them (so \f$\rho\f$
//! and the \f$\Phi\f$ tables are built once per ITERATION for both terms, not once per term) are decided
//! here -- they are implementation, and a builder that has to name them is a builder that can pair them
//! wrong.  \a polarized picks the spin-native pair (tier 4b); \a exch must then be channel-native.
//! \tparam C the concrete correlation functional -- it must satisfy BOTH faces (the unpolarized term takes
//! the plain \c ExFunctional, the polarized one the two-channel \c SpinCorrelation), which is exactly what
//! lets one call serve both branches.
template <class C> std::vector<std::unique_ptr<cDynamic_HT>>
MakeVxcTerms(const std::shared_ptr<ExFunctional>& exch, const std::shared_ptr<C>& corr,
             const std::shared_ptr<const BasisSet::cFIT_SF_ABS>& fb, bool polarized,
             BasisSet::FitQuadrature quad={})
{
    std::shared_ptr<const XC_Quadrature> q=MakeXCQuadrature(fb, std::move(quad));
    std::vector<std::unique_ptr<cDynamic_HT>> terms;
    if (polarized)
    {
        terms.push_back(std::make_unique<Vxc_QuadraturePol  >(exch, q));
        terms.push_back(std::make_unique<Vcorr_QuadraturePol>(corr, q));
    }
    else
    {
        terms.push_back(std::make_unique<Vxc_Quadrature>(exch, q));
        terms.push_back(std::make_unique<Vxc_Quadrature>(corr, q));
    }
    return terms;
}

} //namespace
