// File: XCPolicy.C  The XC-quadrature RUN POLICY: the Becke recipe, the cost model, and the
// Uniform-vs-Becke selector.
//
// Split from qchem.Mesh deliberately.  qchem.Mesh is the geometry-free VALUE layer (Mesh, MeshParams, the
// Nyquist arithmetic) -- no I/O, no environment, no opinions.  Everything here is POLICY: which grid a run
// should use, what the calibrated resolution is, and what to tell the user when they overrule it.  Keeping
// the two apart is what lets the value layer stay quiet and the policy layer stay loud.
//
// WHY THE POLICY IS HERE AT ALL (D6, 2026-08-07): the production Becke recipe and the Auto resolution used
// to live in the anonymous namespace of IntegrationTests/GPW_SCF_UT.C -- a test file defining production
// behaviour.  Tests read this policy now; they do not define it.
//
// THE SELECTOR (V1.26, user ruling 2026-08-07).  Uniform vs Becke is a SMOOTHNESS question, and smoothness
// turns out to be a COST question with a computable crossover, because the two grids scale differently:
//
//   uniform:  n^3, n = UniformDivisions(a, E)          =>  ~ a^3 * alpha_max^(3/2)   -- grows with the
//                                                          cell volume AND the sharpness of the basis
//   Becke:    nAtoms * nRadial * nDirs(scheme, degree) =>  INDEPENDENT of alpha_max  -- grows with nAtoms
//
// So for a soft PP in a big smooth cell uniform wins easily, for sharp Gaussians Becke wins easily, and the
// crossover is a division rather than a judgement call.  Measured with the production recipe: Si/sipp
// (alpha_max=2, a=10.26) wants 6,859 uniform points against ~36,000 free-Becke; F (alpha_max=40, a=8.7)
// wants 357,911 against the same ~36,000.  A factor ~7 either way -- which is why the user's "GPW with any
// high exponents have to, in practice, use Becke" is a statement about cost, not taste.
//
// THREE THINGS THE COST MODEL IS NOT ALLOWED TO FORGET (all one asymmetry -- see kUniformMargin):
//   (1) The comparison is NOT equal-accuracy.  UniformDivisions sizes for the band-limited DENSITY, but the
//       XC integrand v_xc(rho) ~ rho^(1/3) is pointwise-NONLINEAR and not band-limited at all (this is what
//       relCutoff exists for).  So the uniform number is a LOWER BOUND on what uniform really needs, while
//       Becke's radial clustering handles the core natively.  The comparison is biased toward uniform.
//   (2) => therefore a MARGIN, not a tie-break: uniform must win by kUniformMargin to be chosen, so the grey
//       zone defaults to the safe grid.
//   (3) => therefore the diagnostics are ASYMMETRIC: forcing Uniform where Becke is cheaper is usually the
//       dearer AND weaker choice and warns; forcing Becke where uniform is cheaper is usually just wasteful,
//       and only says so.
//
// "USUALLY", NOT "ALWAYS" -- a counterexample the tree already knew about, found by running these diagnostics
// over the suite (GPW_SCF.SiPseudoAtomInBoxMatchesFinite).  BECKE IS NOT UNCONDITIONALLY SAFE: a freely
// rotating DEGENERATE density (that gate's half-filled 3p atom) has orientation-DEPENDENT quadrature error on
// Becke's FIXED-AXIS angular grid -- v_xc is nonlinear, so an anisotropic rho's error rotates with it -- which
// turned an energy-neutral rotation into a ~Ha-scale oscillation.  That test pins Uniform for exactly this
// reason.  So the warning states the COST FACT and names the exception; it does not prescribe Becke.
// (Smearing fixes it -- fractional occupation restores the symmetric density -- so the exception is narrow,
// but it is real, and it is the one case where the asymmetry above inverts.)
module;
#include <cassert>
#include <cstdlib>   // std::getenv/std::atoi/std::atof (the GPW_BECKE_* sweep instruments)
#include <iostream>
#include <string>
export module qchem.Mesh.XCPolicy;
export import qchem.Mesh;
export import qchem.Mesh.Angular;   // LebedevMenu (the degree -> direction-count table)

export namespace qchem::qcMesh
{

//! \brief THE calibrated periodic-Becke XC-quadrature recipe -- the production default for a DFT lattice run.
//!
//! Gate-calibrated on Si (dExc=1.1e-4, dVxc=3.5e-4 against Ecut=60); the angular ladder measured 2026-07-30
//! puts GL-11 at the sub-mHa sweet spot with GL-17/29 on the comparison floor.  Any argument left at its
//! sentinel (\a nRadial<0, \a mhlAlpha<0, \a angularDegree<0) takes the calibrated value, which the matching
//! environment variable may override; an argument passed EXPLICITLY is pinned and ignores the environment,
//! so a probe that has published a resolution keeps it under a sweep.
//!
//! \warning These two defaults carry a lot of weight: they set the ENTIRE Becke side of the selector's cost
//! comparison (nAtoms x nRadial x nDirs), so an over-generous recipe moves the Uniform-vs-Becke crossover
//! against Becke.
//!
//! \par MEASURED (V2.6, 2026-08-07) on THREE systems -- Si diamond (covalent), NaF rocksalt (ionic), Mn
//! sextet atom-in-box (open-shell d).  Instrument: \c GPW_SCF.DISABLED_BeckeRecipeLadder_* (hand-run) -- one
//! converged density, quadratured on a ladder of meshes against a fine reference in the same family, scored
//! on E_xc and the V_xc MATRIX (never \f$\Delta E_{total}\f$ -- the D8 pin).  Verdicts, at the Becke gate's
//! own \c max|dVxc|<1e-3:
//!   - ANGULAR: **degree 29 is ~2.8x over-generous.**  Degree needed: Mn 9, Si 11-15, NaF 15-17.  Recommended
//!     default \b 17 (162 directions vs 450), which leaves NaF -- the worst case -- at 2.2e-4, a 4.5x margin.
//!     Degree 15 is NOT recommended despite sufficing on Si: NaF sits at 7.6e-4, inside tolerance by 1.3x.
//!   - RADIAL: **40 is right, and all three systems agree.**  30 is the first rung inside tolerance
//!     everywhere; 20 is 4-12x out; 25 is marginal on two of the three.
//!   The default is deliberately LEFT AT 29 -- flipping it moves every pinned GPW anchor, so it is a
//!   separate, deliberate re-pin (V2.6a).
//!
//! \par What actually drives the angular requirement (measured, and it is NOT what one would guess):
//! the BECKE PARTITION between DISSIMILAR neighbours, not the on-site density's asphericity.  Mn is a SINGLE
//! atom in a box -- no interatomic partition at all -- and is trivially easy (degree 9) despite being the
//! "most aspherical" system nominally; NaF, two ions of very different size and sharpness, is the HARDEST.
//! Si (partition between IDENTICAL atoms) sits between.  The prediction that ionic would be the most
//! forgiving, being made of near-spherical closed-shell ions, is REFUTED -- it is the least.  The codebase
//! already recorded the mechanism (\c src/Structure/tests/MolecularMeshTests.C: "the fuzzy Voronoi switching
//! shell is angular-quadrature limited"); the site-asphericity framing simply ignored it.
//!
//! \note A cheap a-priori RADIAL check does work, but not the obvious one.  Counting radial nodes inside the
//! sharpest density feature rates \c nRadial=20 as comfortable where measurement puts it 4-50x out -- MHL
//! clusters as \f$r\propto x^m\f$, so those nodes bunch at tiny \a r and leave a gap exactly at the density
//! peak.  The quantity that DOES track the measurement across all three systems is the LOCAL SPACING at the
//! peak of \f$r^2e^{-2\alpha_{\max}r^2}\f$: \f$r_{peak}/\Delta r(r_{peak})\gtrsim3\f$ marks the first rung
//! inside tolerance in every case (Si 3.4 at nR=30, Mn 3.0 at nR=40).  See V2.6 in doc/CleanupCandidates.md.
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

//! \brief What the Uniform-vs-Becke comparison needs to know about a run: its geometry and its two
//! SHARPNESS sources, both expressed as Gaussian exponents so one predicate covers both.
//!
//! A default-constructed value means "no information", which the selector reads as "choose the safe grid"
//! rather than guessing -- see \c ResolveXCMesh.
struct XCMeshSharpness
{
    double cellEdge = 0.0;   //!< Longest cell edge (a.u.).  Binds the isotropic uniform division count.
    int    nAtoms   = 0;     //!< Atoms in the cell.  The Becke mesh is one atom-centred grid per atom.
    //! \brief Sharpest orbital primitive exponent \f$\alpha_{\max}\f$ (\c LatticeSum1E::MaxExponent).  The
    //! DENSITY's sharpest feature is the product of the two tightest primitives, exponent
    //! \f$2\alpha_{\max}\f$ -- which is what \c cutoffFactor already calibrates.
    double alphaMax = 0.0;
    //! \brief Sharpest LOCAL-PP Gaussian exponent \f$\alpha_{pp}=1/2r_{loc}^2\f$
    //! (\c LocalPotential_Gaussian::ShortRangeGaussian).  0 = no PP, or a model with no closed-Gaussian
    //! short part (which cannot be measured -- do not silently treat that as "smooth").  It matters because
    //! the PP integrand \f$\langle\chi_i|V_{short}|\chi_j\rangle\f$ carries exponent
    //! \f$2\alpha_{\max}+\alpha_{pp}\f$, i.e. a sharp PP raises the uniform grid's requirement even under a
    //! soft basis.  Shipped GTH spread: Na q9 15.4, F 10.5, Si 2.6, Na q1 0.64.
    double alphaPP  = 0.0;
    //! \brief Is the run's space group imposed?  The site-adapted invariant Becke mesh prices ~2x the free
    //! mesh (measured 1.97x at degree 29), so it must enter the Becke side of the comparison.
    bool   imposed  = false;

    bool Known() const {return cellEdge>0.0 && nAtoms>0 && alphaMax>0.0;}
};

//! \brief The margin uniform must win by to be selected.  NOT a tuning knob for accuracy -- it exists
//! because the cost comparison is systematically biased toward uniform (see the file header), so a bare
//! `<` would pick the riskier grid throughout the grey zone.
inline constexpr double kUniformMargin = 2.0;

//! \brief The density-scale cutoff a uniform mesh must reach for this run: \f$C\alpha_{\max}+\alpha_{pp}\f$.
//!
//! The first term is the density product \f$2\alpha_{\max}\f$ as \c cutoffFactor already calibrates it; the
//! second is the local PP, which sharpens the \f$\langle\chi|V|\chi\rangle\f$ integrand beyond the density.
double RequiredUniformCutoff(const XCMeshSharpness&, double cutoffFactor=2.0);

//! \brief Point count of the uniform mesh that meets \c RequiredUniformCutoff -- \f$n^3\f$.
long UniformMeshCost(const XCMeshSharpness&, double cutoffFactor=2.0);

//! \brief Point count of the atom-centred Becke mesh \a mp builds: \f$n_{atoms}n_{radial}n_{dirs}\f$.
//!
//! An OVER-estimate by design: the builder drops ~30% of points in the diffuse tails, and over-stating the
//! cost of the safe grid biases the selector toward it, which is the direction we want to err in.  Doubled
//! for an imposed run (the site-adapted invariant mesh).
long BeckeMeshCost(const MeshParams&, const XCMeshSharpness&);

//! \brief Directions a scheme delivers at a requested degree: \f$(L+1)^2/2\f$ for GaussLegendre, the
//! cheapest sufficient tabulated rule for Lebedev.  QUIET -- unlike \c ResolveLebedev it announces nothing,
//! because a cost probe must not narrate a rounding the run may never perform.
int AngularDirections(AngularKind, int degree);

//! \brief Resolve the XC-quadrature choice for a run whose sharpness is KNOWN.
//!
//! - explicit \c Uniform or \c Becke: HONOURED, with an ASYMMETRIC diagnostic -- a warning for the
//!   strictly-dominated choice (forced Uniform where Becke is cheaper, or a Uniform mesh too coarse for the
//!   system), an information line for the merely-wasteful one.  This half is LIVE.
//! - \c cellKind==Auto: the cost comparison is computed and ANNOUNCED, but the answer is still the calibrated
//!   Becke recipe unless the selector is ARMED (\c GPW_XCGRID_AUTOSELECT=1).
//!
//! \warning WHY \c Auto DOES NOT YET ACT ON THE VERDICT.  The margin is not calibrated, and the cost model
//! disagrees with a measurement already in the tree: for Si/sipp it makes uniform ~3-5x cheaper, which would
//! undo the 2026-08-01 Becke-default flip that was itself the result of a measurement.  Since the comparison
//! is known to be biased toward uniform (it sizes for the band-limited density, not for the nonlinear
//! \f$v_{xc}(\rho)\f$ -- see the file header), the disagreement is exactly what \c kUniformMargin is supposed
//! to absorb, and its value is a guess until the D8-compliant grid-convergence measurement fixes it.  Arming
//! it before then would move every pinned anchor on a heuristic.  Announcing it on every run makes each
//! existing GPW test print its own crossover data, which IS the instrument that calibration needs.
MeshParams ResolveXCMesh(const MeshParams& mp, const XCMeshSharpness&);

//! \brief The no-information overload: resolve \c Auto to Becke, the SAFE grid.
//!
//! Not a convenience -- a deliberate fallback for a caller that cannot state the run's sharpness.  Becke is
//! correct for every system and merely wasteful for smooth ones, so it is the right answer to "I do not
//! know"; guessing uniform from a default-constructed sharpness would be an under-resolution.
MeshParams ResolveXCMesh(const MeshParams& mp);

} //export namespace qchem::qcMesh

//-----------------------------------------------------------------------------------------------
namespace qchem::qcMesh
{

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

int AngularDirections(AngularKind kind, int degree)
{
    assert(degree>=0);
    if (kind==AngularKind::GaussLegendre)
        return ((degree+1)/2)*(degree+1);   // n_theta * n_phi, exactly as GaussLegendreAngular builds it
    // Lebedev: the cheapest tabulated rule of at least this degree.  Walks the menu directly rather than
    // calling ResolveLebedev, whose job is to ANNOUNCE the rounding -- correct for a real build, wrong for a
    // cost probe on a grid we may not end up building.
    const LebedevRule* best=nullptr;
    for (const auto& r : LebedevMenu())
        if (r.degree>=degree && (!best || r.nDir<best->nDir)) best=&r;
    return best ? best->nDir : LebedevMenu().back().nDir;   // past the menu: price the densest rule we have
}

double RequiredUniformCutoff(const XCMeshSharpness& s, double cutoffFactor)
{
    assert(s.alphaMax>0.0 && "RequiredUniformCutoff: no basis sharpness to size against");
    return cutoffFactor*s.alphaMax + s.alphaPP;
}

long UniformMeshCost(const XCMeshSharpness& s, double cutoffFactor)
{
    const long n=UniformDivisions(s.cellEdge, RequiredUniformCutoff(s,cutoffFactor));
    return n*n*n;
}

long BeckeMeshCost(const MeshParams& mp, const XCMeshSharpness& s)
{
    assert(s.nAtoms>0 && "BeckeMeshCost: a Becke mesh is one atom-centred grid per atom");
    const long dirs=AngularDirections(mp.angular, mp.angularDegree);
    return long(s.nAtoms)*mp.nRadial*dirs*(s.imposed ? 2 : 1);
}

// The one place the choice is made and narrated.  Both branches print, because the standing pin is that a run
// always says which XC route is in play -- and here it must also say WHY, since the answer now depends on the
// system rather than on a fixed default.
MeshParams ResolveXCMesh(const MeshParams& mp, const XCMeshSharpness& s)
{
    if (!s.Known()) return ResolveXCMesh(mp);   // caller could not state the sharpness -> the safe grid

    const MeshParams becke  =BeckeXCParams();
    const long       cBecke =BeckeMeshCost(becke,s);
    const long       cUnif  =UniformMeshCost(s);
    const double     eReq   =RequiredUniformCutoff(s);
    const bool       uniformWins = double(cUnif)*kUniformMargin <= double(cBecke);

    if (mp.cellKind==UnitCellKind::Auto)
    {
        // ALWAYS announce the comparison: this line is the calibration instrument (see the interface warning).
        std::cout<<"[XC grid choice] Auto: uniform "<<cUnif<<" pts (eCut="<<eReq<<" Ha) vs Becke "<<cBecke
                 <<" pts  [alpha_max="<<s.alphaMax<<", alpha_pp="<<s.alphaPP
                 <<(s.imposed?", imposed":", free")<<"]  => selector would pick "
                 <<(uniformWins?"UNIFORM":"BECKE")<<std::endl;
        static const bool armed=[]{ const char* e=std::getenv("GPW_XCGRID_AUTOSELECT"); return e && std::atoi(e)!=0; }();
        if (!armed || !uniformWins)
        {
            if (armed) std::cout<<"[XC grid choice] Auto -> BECKE (selector armed, Becke wins)"<<std::endl;
            return becke;   // the calibrated default until the margin is measured
        }
        // SIZE the uniform mesh from the basis.  Handing back a bare cellKind would leave nUniform's
        // basis-blind default of 20 in charge, which is the under-resolution this selector exists to prevent.
        MeshParams u=mp;
        u.cellKind=UnitCellKind::Uniform;
        u.eCut    =eReq;
        std::cout<<"[XC grid choice] Auto -> UNIFORM (selector ARMED via GPW_XCGRID_AUTOSELECT; eCut="
                 <<eReq<<" Ha)"<<std::endl;
        return u;
    }

    // EXPLICIT: honoured either way.  The two mistakes are not symmetric, so neither are the two messages --
    // Becke is never WRONG, only sometimes wasteful, while an under-resolved uniform grid is wrong.
    if (mp.cellKind==UnitCellKind::Uniform)
    {
        // Is the mesh the caller actually specified fine enough?  This is the question `nUniform` cannot
        // answer on its own, which is what UniformCutoff is for.
        const double eHave = mp.eCut>0.0 ? mp.eCut : UniformCutoff(s.cellEdge, mp.nUniform);
        if (eHave < eReq)
            std::cerr<<"[XC grid] WARNING: the requested UNIFORM XC mesh resolves eCut="<<eHave
                     <<" Ha but this system needs "<<eReq<<" Ha (alpha_max="<<s.alphaMax
                     <<", alpha_pp="<<s.alphaPP<<"): rho is under-resolved on it -- charge leaks off-grid and "
                       "the XC energy is quadratured on an aliased density.  Raise eCut/nUniform (the remedy "
                       "that is always right), or switch to UnitCellKind::Becke (estimated "<<cBecke
                     <<" pts) IF the system suits it -- see the degenerate-density caveat on the sibling "
                       "warning."<<std::endl;
        if (cBecke < cUnif)
            std::cerr<<"[XC grid] WARNING: UNIFORM was requested but is the DEARER choice here -- resolving "
                       "this system on it WOULD NEED "<<cUnif<<" pts against Becke's "<<cBecke
                     <<", and a uniform grid is also the weaker quadrature for the pointwise-nonlinear "
                       "v_xc(rho).  Becke is usually the better trade here -- but check it suits the system "
                       "first (a freely-rotating DEGENERATE density gets orientation-dependent error on "
                       "Becke's fixed-axis angular grid)."<<std::endl;
        return mp;
    }
    assert(mp.cellKind==UnitCellKind::Becke);
    if (double(cUnif)*kUniformMargin <= double(cBecke))
        std::cout<<"[XC grid] note: BECKE was requested ("<<cBecke<<" pts) where a uniform mesh would cost "
                 <<cUnif<<" -- honoured, and never wrong, just denser than this smooth system needs."<<std::endl;
    return mp;
}

MeshParams ResolveXCMesh(const MeshParams& mp)
{
    if (mp.cellKind!=UnitCellKind::Auto) return mp;
    return BeckeXCParams();   // the calibrated default recipe (nR=40, alpha=2, GL degree 29; GPW_BECKE_* sweepable)
}

} //namespace qchem::qcMesh
