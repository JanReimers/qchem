// File: Mesh.C  Clean-room replacement for qcMesh.  The geometry-free Mesh value type.
//
// A Mesh is a concrete VALUE type = quadrature points + weights, stored SEPARATELY (SoA).
// Points are streamed by evaluators (phi(r)); weights by integrators.  They live in totally
// different algorithms, so they get totally separate arrays.  No polymorphism, no Clone.
//
// Everything lives in namespace qchem::qcMesh so the new library coexists with the old qchem.Mesh
// (which exports a global `class Mesh`) during the migration; the namespace goes away once the
// old library is deleted.
module;
#include <utility>
#include <cassert>
#include <cstdlib>   // std::getenv/std::atoi/std::atof (the GPW_BECKE_* sweep instruments)
#include <string>
export module qchem.Mesh;
export import qchem.Types;

export namespace qchem::qcMesh
{

//! \brief A quadrature mesh: points r_i and weights w_i stored as separate arrays (SoA).
//!
//! \f$\intop f\,d^3r\approx\sum_i w_i f(r_i)\f$.  This is a plain value type with no virtuals:
//! builders \c Append points, consumers read \c Points() (in op(r) evaluators) and \c Weights()
//! (in the integration algorithms) independently.
class Mesh
{
public:
    Mesh() = default;
    //! Build from ready-made SoA arrays (the efficient path when the size is known up front,
    //! e.g. ProductMesh = nRadial*nAngular).  Sizes must match.
    Mesh(rvec3vec_t r, rvec_t w) : itsR(std::move(r)), itsW(std::move(w)) {}

    const rvec3vec_t& Points () const {return itsR;} //!< r_i, for the phi(r) evaluators.
    const rvec_t&     Weights() const {return itsW;} //!< w_i, for the integrators.
    size_t            size   () const {return itsW.size();}

    //! Append one (point,weight) pair.  Builders (the incremental Becke mesh) push here.
    void Append(const rvec3_t& r, double w);
    //! Translate every point: r_i += o.  Used to place a single-center mesh at an atom.
    void ShiftOrigin(const rvec3_t& o);

private:
    rvec3vec_t itsR;
    rvec_t     itsW;
};

//! \brief Radial mesh family.  See the per-class transplanted formulae in Internal/.
enum class RadialKind  {MHL, Log, Linear};
//! \brief Angular mesh family.  All schemes are normalised so \f$\sum_i w_i = 4\pi\f$.
// EulerMaclaren RETIRED 2026-08-07 (user): it was a CONVERGENCE rule, not an exactness rule -- its
// theta quadrature is a transformed trapezoid, so it has NO algebraic degree at all (measured -1: it
// could not integrate even the constant to 1e-10, its weights summing to 4pi only approximately).
// Nothing selected it in production, and GaussLegendre already provides arbitrary L WITH exactness at
// the same ~L^2/2 product-grid cost.  Both survivors now have a real polynomial degree, which is what
// lets the angular knob become degree-typed (R2.15).
enum class AngularKind {Lebedev, GaussLegendre};
//! \brief Which quadrature a periodic UnitCell builds: the uniform FFT-compatible midpoint grid
//! (Hartree/PP -- resolution from \c eCut / \c nUniform), or the atom-centred periodic Becke
//! fuzzy-Voronoi grid (dense radial near each nucleus, cheap diffuse tails -- the XC quadrature;
//! uses the \c radial / \c angularDegree / \c beckeOrder knobs).  See doc/GPWPlan1.md "Becke XC grid".
//! \brief Which quadrature a UnitCell builds for a lattice mesh.  \c Auto = "the caller did not choose",
//! resolved by \c ResolveXCMesh (below) to the calibrated periodic-Becke recipe.  Any consumer reached
//! WITHOUT resolution treats \c Auto as the historical \c Uniform (every consumer compares \c ==Becke, so
//! Auto falls through safely) -- so resolve at the point the mesh spec enters the Hamiltonian.
enum class UnitCellKind {Uniform, Becke, Auto};

//! \brief Becke's iterated smoothing polynomial mapped to the cell cutoff:
//! \f$ s(\mu) = \tfrac12(1-f_k(\mu)),\; f(\mu)=\tfrac12(3\mu-\mu^3) \f$ applied \f$k+1\f$ times
//! (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)).  \f$ s:[-1,1]\to[1,0] \f$, with \f$s(\mp1)\f$
//! pinned at 1/0 and the transition sharpened by each iteration.  Shared by the molecular
//! (MakeMolecularMesh) and periodic (UnitCell) fuzzy-Voronoi partitions.
double BeckeCutoff(double mu, int k);

//! \brief Typed, fully-defaulted mesh parameters.  Set only the knobs you actually use
//! (C++20 designated initializers); no field is required-but-unused.
struct MeshParams
{
    RadialKind  radial    = RadialKind::MHL;    int    nRadial   = 30;
    int         mhl_m     = 2;                   double mhl_alpha = 1.0;     //!< MHL only.
    double      logStart  = 1.0e-4;             double logStop   = 50.0;     //!< Log only.
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
    int         beckeOrder= 3;   //!< Becke fuzzy-Voronoi smoothing iterations (molecular mesh only).
    int         nUniform  = 20;  //!< Uniform periodic real-space grid: points per cell axis (lattice mesh only; \f$n^3\f$ total). Manual fallback when \c eCut<=0.
    double      eCut      = 0.0;  //!< Real-space integration-mesh energy cutoff (a.u.). If >0, the uniform lattice mesh DERIVES its \c nUniform from the Nyquist bound \f$n\gtrsim 2a\sqrt{2E_{cut}}/\pi\f$ (\f$a\f$=longest cell edge; the \f$\times2\f$ is the density bandwidth), and \c nUniform is ignored. 0=use the manual \c nUniform.
    double      relCutoff = 1.0;  //!< Fit-grid density multiplier (CP2K \c REL_CUTOFF): the fit basis scales its \f$E_{cut}\f$ by this. 1=wavefunction bandwidth (LDA); GGA wants >1. Set by the Hamiltonian from the functional's \c GridCutoffFactor().
    UnitCellKind cellKind = UnitCellKind::Uniform;  //!< Lattice mesh only: which quadrature a UnitCell builds (uniform midpoint grid vs periodic Becke).

    //! \brief Compact, deterministic identity string for these parameters.  Two MeshParams give the
    //! same ID() iff they build the same quadrature, so it is the cache key for any mesh-quadrature
    //! quantity that must NOT be shared across different meshes (e.g. a fit basis's numerical Norm --
    //! same basis, different mesh => different normalisation).  All fields are folded in; over-keying
    //! on an inactive field (e.g. logStart when radial=MHL) is harmless.
    std::string ID() const
    {
        using std::to_string;
        return "M{r"  + to_string(static_cast<int>(radial))  + ",nr" + to_string(nRadial)
             + ",m"   + to_string(mhl_m)    + ",a"  + to_string(mhl_alpha)
             + ",ls"  + to_string(logStart) + ",le" + to_string(logStop)
             + ",ang" + to_string(static_cast<int>(angular)) + ",ad" + to_string(angularDegree)
             + ",ar" + to_string(angRot)
             + ",bo" + to_string(beckeOrder)
             + ",nu"  + to_string(nUniform) + ",ec" + to_string(eCut)
             + ",rc" + to_string(relCutoff) + ",ck" + to_string(static_cast<int>(cellKind)) + "}";
    }
};

//! \brief THE calibrated periodic-Becke XC-quadrature recipe -- the production default for a DFT lattice run.
//!
//! Gate-calibrated on Si (dExc=1.1e-4, dVxc=3.5e-4 against Ecut=60); the angular ladder measured 2026-07-30
//! puts GL-11 at the sub-mHa sweet spot with GL-17/29 on the comparison floor.  Any argument left at its
//! sentinel (\a nRadial<0, \a mhlAlpha<0, \a angularDegree<0) takes the calibrated value, which the matching
//! environment variable may override; an argument passed EXPLICITLY is pinned and ignores the environment,
//! so a probe that has published a resolution keeps it under a sweep.
//!
//! Environment instruments (sweep a whole run without rebuilding):
//!   - \c GPW_BECKE_NR     radial point count.                                        Default 40.
//!   - \c GPW_BECKE_ALPHA  MHL radial scale (smaller = nodes pulled toward the core). Default 2.0.
//!   - \c GPW_BECKE_L      angular POLYNOMIAL DEGREE, one meaning for both schemes (R2.15).  Default 29.
//!   - \c GPW_BECKE_ANG    \c "lebedev" selects the Lebedev tables, which resolve the requested degree to
//!                         the cheapest tabulated rule delivering at least it; anything else = GaussLegendre,
//!                         which takes the degree directly.  Because both schemes now read the knob as a
//!                         degree, this A/B is like-for-like (measured on Si: Lebedev beats same-degree GL on
//!                         V_xc elements via the O_h-orbit cancellation, but is 5-10x worse on rho-weighted
//!                         integrals -- its \f$\langle111\rangle\f$ orbit sits on the diamond bond axes).
//!   - \c GPW_BECKE_ROT    radians; rigid generic rotation of the angular grid, which steers those special
//!                         orbits off the bond axes (doc/SymmetryUpgradePlan.md §6a; free runs only).
MeshParams BeckeXCParams(int nRadial=-1, double mhlAlpha=-1.0, int angularDegree=-1);

//! \brief Resolve \c UnitCellKind::Auto to the actual XC quadrature; an explicit \c cellKind is always honored.
//!
//! Auto means "the caller did not choose", and the answer is the atom-centred periodic BECKE mesh
//! (\c BeckeXCParams above) -- the near-ideal grid for the one pointwise-nonlinear sharp-at-core term,
//! diffuse bases included (the 2026-08-01 Becke-default flip, doc/GPWPlan1.md).
//!
//! This resolution consults NO run context, which is what lets it live here beside the enum rather than in a
//! driver: the last context-dependent branch was the BZ-reduced carve-out, retired 2026-08-02 (plan §7 step 5)
//! once the mixed-rule site-adapted INVARIANT Becke mesh was gate-verified -- it prices ~2x the free mesh at
//! the production recipe (measured 1.97x at degree 29), so imposed runs get Becke too instead of falling back
//! to the uniform raster.  Should a future policy genuinely need the run context, it belongs in a
//! \c SymmetryPolicy-style object and this function becomes its default.
MeshParams ResolveXCMesh(const MeshParams& mp);

} //export namespace qchem::qcMesh

//-----------------------------------------------------------------------------------------------
namespace qchem::qcMesh
{

// Append grows the SoA arrays by one.  blaze resize copies, so this is a setup-path convenience
// (the incremental Becke builder), NOT for hot loops -- when the size is known, use the
// from-arrays constructor instead.
void Mesh::Append(const rvec3_t& r, double w)
{
    size_t n=itsW.size();
    assert(itsR.size()==n);
    itsR.resize(n+1,true);
    itsW.resize(n+1,true);
    itsR[n]=r;
    itsW[n]=w;
}

void Mesh::ShiftOrigin(const rvec3_t& o)
{
    for (size_t i=0; i<itsR.size(); i++) itsR[i]+=o;
}

double BeckeCutoff(double mu, int k)
{
    for (int i=k; i>=0; i--) mu=0.5*(3*mu - mu*mu*mu);
    return 0.5*(1.0-mu);
}

MeshParams BeckeXCParams(int nRadial, double mhlAlpha, int angularDegree)
{
    auto envi=[](const char* n, int    d){ const char* s=std::getenv(n); return s ? std::atoi(s) : d; };
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };
    if (nRadial      <0)   nRadial      =envi("GPW_BECKE_NR",    40);
    if (mhlAlpha     <0.0) mhlAlpha     =envd("GPW_BECKE_ALPHA", 2.0);
    if (angularDegree<0)   angularDegree=envi("GPW_BECKE_L",     29);
    MeshParams mp;
    mp.cellKind=UnitCellKind::Becke;
    mp.radial =RadialKind::MHL;            mp.nRadial =nRadial; mp.mhl_m=2; mp.mhl_alpha=mhlAlpha;
    const char* ang=std::getenv("GPW_BECKE_ANG");
    mp.angular=(ang && std::string(ang)=="lebedev") ? AngularKind::Lebedev : AngularKind::GaussLegendre;
    mp.angularDegree=angularDegree;   // ONE meaning for both schemes: GL takes it directly, Lebedev resolves
                                      // it to the cheapest rule of at least that degree (R2.15).
    mp.angRot=envd("GPW_BECKE_ROT", 0.0);
    return mp;
}

MeshParams ResolveXCMesh(const MeshParams& mp)
{
    if (mp.cellKind!=UnitCellKind::Auto) return mp;
    return BeckeXCParams();   // the calibrated default recipe (nR=40, alpha=2, GL degree 29; GPW_BECKE_* sweepable)
}

} //namespace qchem::qcMesh
