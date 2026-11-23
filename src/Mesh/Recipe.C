// File: Mesh/Recipe.C  qchem.Mesh.Recipe -- what a mesh is TO BE (MeshParams and its kinds), apart from the BUILT point set (qchem.Mesh).
//
// V-MESH-FACTORY step 1 (user ruling 2026-11-21, "Please proceed with the plan"; lane H design note §3): a RECIPE churns -- a new
// radial kind, a partition kind, a Becke knob -- while the built Mesh every hub holds rarely changes.  Sharing one module made every
// recipe edit rebuild the ~500 transitive importers of qchem.Mesh.  Step 1 split it out; step 5 ended qchem.Mesh's re-export, so
// every importer states it (507 -> 214 transitive importers of a recipe edit, doc/MeshFactoryPlan.md §6).
module;
#include <string>
export module qchem.Mesh.Recipe;
export import qchem.Types;
export import qchem.PointSampling;   // SamplingAccuracy (MeshParams::sampling)
import qchem.Symmetry.Group;         // Symmetry::GroupPtr / Frame: what a stated point set FORCES (pin 30)

export namespace qchem::qcMesh
{

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
