// File: Mesh.C  Clean-room replacement for qcMesh.  The geometry-free Mesh value type.
//
// A Mesh is a concrete VALUE type = quadrature points + weights, stored SEPARATELY (SoA).
// Points are streamed by evaluators (phi(r)); weights by integrators.  They live in totally
// different algorithms, so they get totally separate arrays.  No Clone, and no subclasses.
// ONE abstract face since 2026-11-04 (user ruling): a Mesh IS a qchem::PointSampling -- the point set a
// Fourier-series field asks HOW to be sampled on it (a NUFFT by default; an exact inverse FFT on the uniform
// cell raster).  The uniform override is a HELD STRATEGY the builder installs, not a subclass, because a Mesh
// is returned and copied by value everywhere and a subclass would slice; so it stays a value type -- a copy
// keeps its sampling.
//
// Everything lives in namespace qchem::qcMesh so the new library coexists with the old qchem.Mesh
// (which exports a global `class Mesh`) during the migration; the namespace goes away once the
// old library is deleted.
module;
#include <utility>
#include <cassert>
#include <string>
#include <vector>    // the site-block index table (Mesh::itsSiteStart)
#include <memory>    // the held FourierSamplingStrategy
export module qchem.Mesh;
export import qchem.Types;
export import qchem.PointSampling;   // the face a Mesh implements (MakeFourierSampler)
import qchem.Symmetry.Group;        // Symmetry::GroupPtr: what a point set FORCES on the field it samples (Group face step 3b)

export namespace qchem::qcMesh
{

class Mesh;
//! \brief THE GEOMETRY DERIVATIVE OF A MESH'S WEIGHTS (lane B 2026-11-19; the XC term's grid response): installed on a mesh by the builder
//! that knows how its weights were made (the periodic Becke build: the cell, the sites, the partition's order and eps), so a consumer asks
//! the mesh, never the structure.  For a per-point field \a f HELD per point index, \c Pullback returns \f$\partial(\sum_g f_gw_g)/\partial R_A\f$
//! per site (the points riding their owners, \c Mesh::OwnerOffsets) and its strain twin \f$\sum_X(\partial/\partial X)\otimes X\f$ over every
//! unfolded position the weights depend on.  A mesh without one has no weight response to give.
class MeshWeightResponse
{
public:
    virtual ~MeshWeightResponse()=default;
    struct Result { rvec3vec_t gradient; rmat3d_t virial; };
    virtual Result Pullback(const Mesh& m, const rvec_t& f) const=0;
};

//! \brief A quadrature mesh: points r_i and weights w_i stored as separate arrays (SoA).
//!
//! \f$\intop f\,d^3r\approx\sum_i w_i f(r_i)\f$.  A value type: builders \c Append points, consumers read
//! \c Points() (in op(r) evaluators) and \c Weights() (in the integration algorithms) independently.  Its one
//! virtual face is \c PointSampling (see the file header): how a Fourier series is sampled on these points.
class Mesh : public virtual PointSampling
{
public:
    Mesh() = default;
    //! Build from ready-made SoA arrays (the efficient path when the size is known up front,
    //! e.g. ProductMesh = nRadial*nAngular).  Sizes must match.
    Mesh(rvec3vec_t r, rvec_t w);
    //! The atom-centred sibling: \a siteStart is the first point index of each SITE BLOCK (see below), and \a ownerOffset (empty, or one
    //! per point) each point's offset from its site's own centre (see \c OwnerOffsets).
    Mesh(rvec3vec_t r, rvec_t w, std::vector<size_t> siteStart, rvec3vec_t ownerOffset={});

    const rvec3vec_t& Points () const override; //!< r_i, for the phi(r) evaluators.
    //! \copydoc PointSampling::MakeFourierSampler
    //! The builder's installed strategy when there is one (the uniform raster's exact inverse FFT), else the NUFFT default.
    std::unique_ptr<FourierSampler> MakeFourierSampler(const FourierBox& b) const override;
    //! A BUILDER records a structure of its points the series can exploit.  \c Append and \c ShiftOrigin drop it (the
    //! points are then no longer the ones it described).
    void InstallFourierSampling(std::shared_ptr<const FourierSamplingStrategy> s);
    //! The builder's WEIGHT RESPONSE (\c MeshWeightResponse), or null when none was installed.  \c Append and \c ShiftOrigin drop it, as
    //! they drop the sampling strategy (the weights are then no longer the ones it differentiates).
    void InstallWeightResponse(std::shared_ptr<const MeshWeightResponse> r);
    const MeshWeightResponse* WeightResponse() const {return itsWeightResponse.get();}
    //! \brief WHAT THESE POINTS FORCE on a field sampled on them, in the Cartesian frame (doc/Pins.md pin 30; Group face step 3b).
    //! Only the builder knows: a uniform cell raster forces its cell's lattice T (a field on it is periodic by construction), the radial-only
    //! set forces O(3) about its one centre (one direction at \f$4\pi\f$: exact only for a spherical integrand), and everything else forces
    //! nothing.  So the builder of a symmetric point set INSTALLS its group; an unstamped mesh is a plain point cloud and answers {e}.
    //! \c Append and \c ShiftOrigin drop the stamp, as they drop the sampling strategy: the points are no longer the ones it described.
    void               InstallForces(Symmetry::GroupPtr g);
    Symmetry::GroupPtr Forces() const;
    //! \copydoc PointSampling::FourierTarget
    //! Every mesh has one, SHARED by its copies (a copy is the same points, so the same consumer's tolerance and calibration).
    SamplingTarget* FourierTarget() const override;
    //! INLINE BY DESIGN (D-BUILD-BLAST, lane B 2026-11-18): \c Weights, \c size and the site-block bounds are read inside point loops,
    //! and the build has no LTO, so an out-of-line body would cost a call per loop test in every importer.  Bodies that a loop does
    //! not call live in Imp/Mesh.C.
    const rvec_t&     Weights() const {return itsW;} //!< w_i, for the integrators.
    size_t            size   () const {return itsW.size();}

    //! \name SITE BLOCKS — the partition an atom-centred mesh already carries
    //! An atom-centred quadrature (Becke and friends) is a UNION of per-centre grids, and its weights
    //! already fold in that centre's partition function \f$w_A(r)\f$: that is what makes
    //! \f$\int f=\sum_A\int w_A f\f$ the scheme's defining identity.  So the mesh knows which centre
    //! generated each point — it just used to throw the fact away, which is why a "site moment" had to be
    //! faked by sampling a field at a guessed point (doc/OpenWork.md Step 0a).  A builder that generates
    //! per centre calls \c MeshBuilder::BeginSite, and then \f$\int w_A f\f$ is EXACTLY the sum of
    //! \f$w_q f_q\f$ over site \f$A\f$'s block — no extra quadrature, no extra parameter, no partition
    //! choice at the consumer.
    //! Blocks are contiguous and cover the mesh; EMPTY (\c NSites()==0) means the mesh has no site
    //! structure (a plain product or uniform grid), and site integrals are then simply unavailable —
    //! callers must ask, not assume.  These carry NO structural types (plain indices), so \c qcMesh stays
    //! a leaf: WHICH atom site \f$A\f$ is belongs to whoever built the mesh.
    //!@{
    size_t NSites  ()          const {return itsSiteStart.size();}
    size_t SiteBegin(size_t a) const {assert(a<NSites()); return itsSiteStart[a];}
    size_t SiteEnd  (size_t a) const {assert(a<NSites()); return a+1<NSites() ? itsSiteStart[a+1] : size();}
    //! \brief Each point's offset \f$u_i\f$ from its OWNER's centre: the point was generated at \f$R_A+u_i\f$ by site \f$A\f$'s template
    //! (before any wrap into a home cell).  What a GEOMETRY DERIVATIVE of an atom-centred quadrature needs: the point rides its owner, so
    //! the partition is re-evaluated at \f$R_A+u_i\f$ and the virial carries \f$u_i\f$ (lane B 2026-11-19).  EMPTY when the builder did not
    //! record it (a product or uniform grid has no owners); otherwise one per point.
    const rvec3vec_t& OwnerOffsets() const {return itsU;}
    //!@}

    //! Append one (point,weight) pair.  Builders (the incremental Becke mesh) push here.
    void Append(const rvec3_t& r, double w);
    //! Translate every point: r_i += o.  Used to place a single-center mesh at an atom.
    void ShiftOrigin(const rvec3_t& o);

private:
    rvec3vec_t          itsR;
    rvec_t              itsW;
    std::vector<size_t> itsSiteStart;   //!< first point index of each site block ({} = no site structure)
    rvec3vec_t          itsU;           //!< each point's offset from its owner's centre ({} = not recorded)
    std::shared_ptr<const FourierSamplingStrategy> itsSampling;   //!< null = the NUFFT default
    std::shared_ptr<const MeshWeightResponse> itsWeightResponse;  //!< null = no weight derivative installed
    Symmetry::GroupPtr  itsForces;      //!< what the builder stamped (null = a plain point cloud: {e})
    std::shared_ptr<SamplingTarget> itsTarget=std::make_shared<SamplingTarget>();   //!< the sampling slot (unset = step-1 tolerance)
};

//! \brief The per-site integrals \f$\int w_A(r)f(r)\,d^3r\f$ of a field \a f sampled on \a m's points —
//! one entry per site block, in the builder's site order.
//!
//! This is the honest way to ask an atom-centred mesh for an ATOMIC quantity (moment, charge, …): it is a
//! genuine integral with units, not a point sample, and it needs no radius and no extra convention beyond
//! the partition the mesh was already built with.  \a f must be sampled on \c m.Points().
//! CAVEAT worth carrying to any consumer: the answer is only as canonical as the partition — a Becke fuzzy
//! basin is a CHOICE.  The partition-free definition is R. F. W. Bader's QTAIM zero-flux basin
//! (\f$\nabla\rho\cdot n=0\f$), a wanted future feature (doc/OpenWork.md Step 0a); until it lands, report
//! WHICH partition produced the number.
rvec_t SiteIntegrals(const Mesh& m, const rvec_t& f);

//! \brief \f$\int f\,d^3r=\sum_g w_g f_g\f$ for a field already SAMPLED on \a m's points -- the tabulated
//! sibling of \c qcMesh::Integrate(const Mesh&, const ScalarFunction<T>&) in \c qchem.Mesh.Quadrature.
//!
//! Lives HERE, beside \c SiteIntegrals, because it is the same shape (a weighted sum over a field the
//! caller already has) and because it must be reachable from anyone holding only a \c Mesh: it is what a
//! consumer of the \c BasisSet::Quadrature face integrates with, whether the mesh behind that face is an
//! FFT raster, a uniform cell grid, or an atom-centred Becke build.
double Integrate(const Mesh& m, const rvec_t& f);

//! \brief Radial mesh family.  See the per-scheme formulae in Internal/.  \c DE2 = Mitani's double-exponential rule (the radial rule of
//! Dasgupta & Herbert's SG-2/SG-3), Internal/DE2RadialMesh.C.
enum class RadialKind  {MHL, Log, Linear, DE2};
//! \brief Angular mesh family.  All schemes are normalised so \f$\sum_i w_i = 4\pi\f$.
// EulerMaclaren RETIRED 2026-08-07 (user): it was a CONVERGENCE rule, not an exactness rule -- its
// theta quadrature is a transformed trapezoid, so it has NO algebraic degree at all (measured -1: it
// could not integrate even the constant to 1e-10, its weights summing to 4pi only approximately).
// Nothing selected it in production, and GaussLegendre already provides arbitrary L WITH exactness at
// the same ~L^2/2 product-grid cost.  Both survivors now have a real polynomial degree, which is what
// lets the angular knob become degree-typed (R2.15).
//! \c None (V-XC-ONE-TERM step 5, 2026-11-14): NO angular rule -- the RADIAL-ONLY point set, \f$\int f(r)\,4\pi r^2dr\f$ about one centre
//! (one direction at weight \f$4\pi\f$).  Exact only for a SPHERICAL integrand, so it imposes O(3) about its centre (H_points = O(3),
//! doc/Pins.md pin 30): admitted only on a structure with full rotation symmetry (the atom path; \c Hamiltonian::AdmitVxc refuses it
//! elsewhere).  \c angularDegree is not read for it (its degree is 0).  It used to be spelled Lebedev at degree 0 (the 1-direction rule),
//! which on the molecule doubles as the UNSTATED sentinel -- so the radial-only set could only be inferred, never stated.
enum class AngularKind {Lebedev, GaussLegendre, None};
//! The polynomial degree an angular rule is EXACT to: 0 for the radial-only set, else the stated \c angularDegree (Lebedev rounds UP).
constexpr int ExactAngularDegree(AngularKind k, int angularDegree) {return k==AngularKind::None ? 0 : angularDegree;}
//! \brief Which quadrature a periodic UnitCell builds: the uniform FFT-compatible midpoint grid
//! (Hartree/PP -- resolution from \c eCut / \c nUniform), or the atom-centred periodic Becke
//! fuzzy-Voronoi grid (dense radial near each nucleus, cheap diffuse tails -- the XC quadrature;
//! uses the \c radial / \c angularDegree / \c beckeOrder knobs).  See doc/GPWPlan1.md "Becke XC grid".
//! \brief Which quadrature a UnitCell builds for a lattice mesh.  \c Auto = "the caller did not choose"; the
//! POLICY layer resolves it (\c qcMesh::ResolveXCMesh in \c qchem.Mesh.XCPolicy, which weighs the two grids'
//! costs against the run's basis and PP sharpness).  Any consumer reached WITHOUT resolution treats \c Auto as
//! the historical \c Uniform (every consumer compares \c ==Becke, so Auto falls through safely) -- so resolve
//! at the point the mesh spec enters the Hamiltonian.
enum class UnitCellKind {Uniform, Becke, Auto};

//! \brief \c MeshParams::cellKind: a \c UnitCellKind that cannot be set to \c Becke by itself (OOD BM(1)).
//!
//! The struct's other defaults are the 30 / 1 / degree-5 grid, not the Becke recipe, so a bare \c MeshParams{.cellKind=Becke} built
//! a degree-5 mesh; that cost a bogus 40 mHa "symmetry bug".  So stating \c Becke does not compile: ask for the recipe,
//! \c BeckeXCParams(), or, where the struct already holds one (the facade's \c xcMesh does) or the caller states every field, say so
//! with \c RecipeStated.  A value known only at run time takes \c RecipeStated too.
class CellKind
{
public:
    constexpr CellKind() = default;
    consteval CellKind(UnitCellKind k) : itsKind(k)
    {
        if (k==UnitCellKind::Becke) throw "cellKind Becke alone: use BeckeXCParams() or CellKind::RecipeStated";
    }
    static constexpr CellKind RecipeStated(UnitCellKind k) {CellKind c; c.itsKind=k; return c;}
    constexpr operator UnitCellKind() const {return itsKind;}
private:
    UnitCellKind itsKind=UnitCellKind::Auto;
};


//! \brief Becke's iterated smoothing polynomial mapped to the cell cutoff:
//! \f$ s(\mu) = \tfrac12(1-f_k(\mu)),\; f(\mu)=\tfrac12(3\mu-\mu^3) \f$ applied \f$k+1\f$ times
//! (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)).  \f$ s:[-1,1]\to[1,0] \f$, with \f$s(\mp1)\f$
//! pinned at 1/0 and the transition sharpened by each iteration.  Shared by the molecular
//! (MakeMolecularMesh) and periodic (UnitCell) fuzzy-Voronoi partitions.
//! INLINE (lane B, 2026-11-08): a derivative leaf.  Enzyme differentiates only bodies it can see, and an importer sees an exported
//! function's body only when it is inline (doc/GeometryOptimizationPlan.md, the AD verdict), so the Becke partition's unit can
//! differentiate through it in place.  The same operations in the same order: bit-identical.
inline double BeckeCutoff(double mu, int k)
{
    for (int i=k; i>=0; i--) mu=0.5*(3*mu - mu*mu*mu);
    return 0.5*(1.0-mu);
}

//! \brief EMPIRICAL: the DE2 radial rule's mapping parameter alpha (\c MeshParams::de2_alpha), one value for every element.
//!
//! There is no rigorous alpha: it only shares a fixed node count between core and tail, and the ε radial climb picks the count.  MEASURED
//! 2026-11-21 (lane A; `~/Code/qchem6-runs/radial_de2/fine/`, `scfrun --radial DE2`): free LDA atoms on the atom route (Slater Medium, the
//! XC radial count stated, Hartree at 800), the fewest points from which E stays within ε of the 300-point value, on a 5-point ladder:
//!   ε 1e-6, alpha 0.6 / 0.8 / 1.0 / 1.2 / 1.5 | MHL (m 3, alpha 2: the atom default)
//!     all-electron  F 50/45/40/45/40 | 40    Na 55/45/50/45/45 | 45    Ar 65/55/55/55/60 | 60    Mn 70/65/65/55/60 | 65
//!     GTH           Na q1 25/30/30/35/35 | 45   Na q9 50/50/45/50/55 | 55   O q6 40/35/35/40/40 | 50   F q7 40/40/40/45/45 | 50
//!                   Si q4 35/25/40/35/45 | 55   Mn q15 45/45/45/50/50 | 55
//!   ε 1e-8 at alpha 1.0 vs MHL: all-electron F 60|60, Na 55|70, Ar 70|75, Mn 80|80; GTH Na q1 40|55, Na q9 60|70, O q6 50|60, F q7 50|60,
//!   Si q4 45|65, Mn q15 55|75.
//! The optimum is FLAT from 0.8 to 1.2 on every atom (within one 5-point rung of the best), so one value serves all elements and a
//! per-element table would hold the same number (user 2026-11-21: "a ball park number to keep grid sizes reasonably down").  1.0 is also
//! Mitani 2011's molecular value; his atom-averaged optimum 1.9 (2012) costs 5-20 points on the GTH atoms here.  A crystal revisits it
//! only if its radial climb shows alpha matters.  DE2 against MHL: equal on all-electron atoms, 15-30% fewer points on pseudo-atoms.
inline constexpr double kDE2AlphaCalibrated = 1.0;

//! \brief Typed, fully-defaulted mesh parameters.  Set only the knobs you actually use
//! (C++20 designated initializers); no field is required-but-unused.
struct MeshParams
{
    RadialKind  radial    = RadialKind::MHL;    int    nRadial   = 30;
    int         mhl_m     = 2;                   double mhl_alpha = 1.0;     //!< MHL only.
    double      logStart  = 1.0e-4;             double logStop   = 50.0;     //!< Log only.
    //! DE2 only: the mapping parameter alpha (a SHAPE parameter: how the nodes divide between core and tail; nRadial sets the accuracy)
    //! and the node span [rMin, rMax] (Mitani's truncation, both ends included).  alpha is a calibration, never part of the formula
    //! (Internal/DE2RadialMesh.C): its default is \c kDE2AlphaCalibrated, whose record is at its declaration.
    double      de2_alpha = kDE2AlphaCalibrated;   double de2_rMin = 1.0e-7;   double de2_rMax = 40.0;
    AngularKind angular   = AngularKind::Lebedev;
    //! \brief The requested angular POLYNOMIAL DEGREE -- the same quantity for every scheme (R2.15).
    //!
    //! It was \c nAngular and meant different things per scheme: a DIRECTION COUNT for Lebedev but the
    //! degree L for the others, so the two could not be compared and the Becke default could not be moved
    //! between them.  (A third meaning, EulerMaclaren's resolution knob, was retired with that scheme --
    //! it had no degree at all.)  Now: Lebedev RESOLVES the degree to the cheapest tabulated rule that
    //! delivers at least it (see \c ResolveLebedev, which announces any rounding); GaussLegendre uses it
    //! directly, as it always did.
    int         angularDegree = 5;      //!< default: the 12-direction Lebedev rule (unchanged grid).
    double      angRot    = 0.0;  //!< Rigid rotation of the angular grid (radians, about the fixed generic axis \f$(1,2,3)/\sqrt{14}\f$).  Quadrature exactness is rotation-invariant; the knob steers a grid's special orbits OFF structure axes (e.g. Lebedev's \f$\langle111\rangle\f$ orbit off diamond's bonds -- doc/SymmetryUpgradePlan.md §6a rotation insight, free runs only).  0 = off (bit-identical historical grids).
    int         beckeOrder= 3;   //!< Becke fuzzy-Voronoi smoothing order, molecular and periodic meshes alike (\c BeckeCutoff).  It iterates order+1 times, so order 2 is Becke's published three iterations.
    int         nUniform  = 20;  //!< Uniform periodic real-space grid: points per cell axis (lattice mesh only; \f$n^3\f$ total). Manual fallback when \c eCut<=0.
    double      eCut      = 0.0;  //!< Real-space integration-mesh energy cutoff (a.u.). If >0, the uniform lattice mesh DERIVES its \c nUniform from the Nyquist bound \f$n\gtrsim 2a\sqrt{2E_{cut}}/\pi\f$ (\f$a\f$=longest cell edge; the \f$\times2\f$ is the density bandwidth), and \c nUniform is ignored. 0=use the manual \c nUniform.
    double      relCutoff = 1.0;  //!< Fit-grid density multiplier (CP2K \c REL_CUTOFF): the fit basis scales its \f$E_{cut}\f$ by this. 1=wavefunction bandwidth (LDA); GGA wants >1. Set by the Hamiltonian from the functional's \c GridCutoffFactor().
    //! Which quadrature the STRUCTURE builds: a lattice resolves \c Auto by policy (\c ResolveXCMesh; unresolved, it reads as
    //! \c Uniform); a molecule or an atom builds its atom-centred mesh for \c Auto or \c Becke and REFUSES a stated \c Uniform
    //! (a cell mesh imposes a translation group the structure lacks).  Default \c Auto (2026-10-31; was \c Uniform, which every
    //! molecular run carried while building Becke anyway -- the default could not be told from a request).
    CellKind cellKind;
    //! \brief Becke competitor-series TOLERANCE (periodic Becke only): the magnitude screen that fixes how many cell images each
    //! partition weight sums.  1e-6 holds the equivalent-site-share gate (site shares equal to 1e-8 relative) with a decade of
    //! margin; 1e-5 breaks it.  A TOLERANCE trade, not a bit-identical restructuring -- the weights move at ~1e-6 relative.
    //! Was the env override GPW_BECKE_EPS read inside the mesh build; typed here (D-ENV step 5).
    double      beckeEps  = 1.0e-6;
    //! \brief The accuracy a Fourier series sampled on this mesh must meet (the rho-tilde sampler's step 2; \c qchem::SamplingAccuracy):
    //! stated by the run for its XC quadrature, unset everywhere else.  It sizes a sampler ON the points, never the points.
    SamplingAccuracy sampling;

};

//-- What a point set FORCES (doc/Pins.md pin 30) -----------------------------------------------------------------------------------
//! \brief WHAT A STATED XC POINT SET FORCES on a field sampled on it, as a Group in \a frame (Group face step 3b; it replaced the three-state
//! \c ImposedSymmetry and \c ImposedBy(MeshParams)).  The KIND fixes the group's type and the frame supplies where; no point is built (user
//! 2026-11-21: "We know any radial forces O(3) without creating lists of r,w pairs"):
//!  - a uniform cell mesh forces T of the frame's lattice (a field on the raster is periodic by construction; it reads no angular rule);
//!  - the radial-only set (\c AngularKind::None: one direction at \f$4\pi\f$, exact only for a spherical integrand) forces O(3) about the
//!    frame's centre;
//!  - an atom-centred (Becke) mesh forces nothing, and so does \c Auto: each structure resolves it to a mesh it admits.
//! It FAILS, with the reason, when the frame lacks what the kind needs (no lattice for a uniform mesh, no centre for the radial-only set):
//! that point set cannot exist there.  ONE source of truth: the structures' mesh builders stamp \c Mesh::Forces from this, in their own
//! frame, and the pin-30 pre-flight (\c Hamiltonian::AdmitVxc) asks it of the structure's group's frame.  Spin symmetry is a SEPARATE axis
//! and never enters (user, 2026-11-16: the real-space group "should not be conflated with spin space").
Outcome<Symmetry::GroupPtr,std::string> Forces(const MeshParams& mp, const Symmetry::Frame& frame);

//-- The uniform-mesh resolution bridge --------------------------------------------------------------
// eCut and nUniform are two spellings of ONE quantity: how finely a uniform midpoint mesh samples the cell.
// The mapping between them was a literal inside UnitCell::CreateIntegrationMesh, so any consumer that wanted
// to ASK the question ("is this mesh fine enough?", "which grid is cheaper?") had to duplicate it -- the
// R2.12 "name it once" case.  It is named here, beside the two fields it relates.

//! \brief Divisions per axis a uniform midpoint mesh over a cell of longest edge \a cellEdge needs in order
//! to resolve a field of cutoff \a eCut.
//!
//! The midpoint rule integrates \f$e^{iG\cdot r}\f$ exactly only while the fractional spacing \f$a/n\f$
//! out-resolves the field's highest \f$G=\sqrt{2E}\f$: \f$n>aG/\pi\f$, doubled for a density-bandwidth
//! field, giving \f$n=\lceil 2a\sqrt{2E}/\pi\rceil\f$.  The longest cell edge binds an isotropic \a n, so a
//! shorter axis is (harmlessly) over-resolved.
int UniformDivisions(double cellEdge, double eCut);

//! \brief The inverse: the cutoff an \a n-per-axis uniform mesh over edge \a cellEdge actually resolves.
//!
//! This is what makes a bare \c nUniform COMPARABLE to an \c eCut -- i.e. what lets a caller ask whether a
//! hand-set division count is fine enough for a given basis, without re-deriving the Nyquist relation.
//! Exact inverse up to \c UniformDivisions' ceiling, so \c UniformCutoff(a,UniformDivisions(a,E)) >= E.
double UniformCutoff(double cellEdge, int n);

} //export namespace qchem::qcMesh
