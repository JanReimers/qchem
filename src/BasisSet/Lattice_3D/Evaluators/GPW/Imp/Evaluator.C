// File: BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C  GPW_Evaluator implementation.
module;
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>    // std::chrono (timing the MakeLocalPP integrate-back at different kappa)
#include <cstdlib>   // std::getenv/std::atof (GPW_LOCALPP_RELCUTOFF / GPW_OMP_THREADS knobs)
module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.Blaze;       // rvec_t, rmat_t, rsmat_t, blazem::zeroH<dcmplx>
import qchem.Vector3D;    // vec3_t + rvec3_t / rvec3vec_t arithmetic (r - R, componentwise add)
import qchem.Mesh.Quadrature; // qcMesh::WeightedOverlap / Overlap (the real-space PP quadrature primitives)
import qchem.ScalarFunction;  // ScalarFunction<double> (the V_loc / beta*Ylm fields handed to the quadrature)
import qchem.Math;            // norm, Pi, sqrt (the real spherical harmonics for the KB projectors)
import qchem.Math.Angular;    // Monomial/CartTerm/SphericalShell (the analytic KB Cartesian expansion)
import qchem.Structure;       // Structure / Atom (the PP centres + Z, and CreateIntegrationMesh)
import qchem.UnitCell;        // UnitCell (the direct cell for CollocateDensity / IntegratePotential grid<->cell)
import qchem.Reporting;       // the run report -- EmitGridsReport builds the `grids` section
import qchem.Symmetry.Lattice_3D.Fold;   // FoldGVectors (the grid report's {G}-star column)

namespace qchem::BasisSet::Lattice_3D
{

namespace
{
// TRIM (k time-reversal-invariant, 2k a reciprocal lattice vector) is no longer re-derived here from the
// float k: the ctor's kIsReal FACT (itsTRIM) comes from the Bloch irrep's Symmetry::IsReal() -- exact
// integer arithmetic on the BZ grid, no tolerance (doc/RealComplexPlan.md Step 1).
//! \brief The Bloch phase \f$e^{2\pi i\,k\cdot n}\f$ of integer cell \a n -- EXACTLY \f$\pm1\f$ at a
//! TRIM \a kFrac.
//!
//! There \f$2k\f$ is integral, so \f$e^{2\pi ik\cdot n}=(-1)^{(2k)\cdot n}\f$: a parity, not a
//! transcendental.  Computing it as \c std::exp leaves \f$\sin(\pm\pi)\approx1.2\times10^{-16}\f$ in
//! the imaginary part, and that junk propagates into S, T, V_loc, the KB projectors and the XC-mesh
//! basis tables at EVERY zone-boundary k -- where it is then symmetrised away, having cost a complex
//! multiply per term.  Exact \f$\pm1\f$ makes those blocks exactly real instead of nearly real, which
//! is the difference between a REAL-ARITHMETIC fast path selectable by a fact and one selectable only
//! by a tolerance somebody has to defend (doc/GPWPlan1.md).  Away from TRIM this is the plain \c exp.
dcmplx BlochPhase(const rvec3_t& kFrac, const ivec3_t& n, bool trim)
{
    if (trim)
    {
        const long p = std::lround(2.0*kFrac.x)*n.x + std::lround(2.0*kFrac.y)*n.y
                     + std::lround(2.0*kFrac.z)*n.z;
        return dcmplx((p & 1L) ? -1.0 : 1.0, 0.0);
    }
    const double kn=kFrac.x*n.x + kFrac.y*n.y + kFrac.z*n.z;
    return std::exp(dcmplx(0.0, 2.0*Pi*kn));
}

// Build a lattice-translation set {R} (Cartesian, origin first) and its matching Bloch phases {e^{ik.R}} for
// a cutoff radius -- the {R}+{phase} weighted point set (future: one qcMesh cMesh).  phase = exp(2 pi i k.n)
// with n the integer cell index (convention-safe).  Rcut<=0 -> the home cell only (origin, phase 1).
void BuildImages(const UnitCell& cell, double Rcut, const rvec3_t& kFrac, bool trim,
                 std::vector<rvec3_t>& R, cvec_t& phase)
{
    blazem::VecBuilder<dcmplx> ph;
    if (Rcut>0.0)
        for (const auto& n : cell.CellsInSphere(Rcut))
        {
            R.push_back(cell.ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z))));
            ph.Append(BlochPhase(kFrac,n,trim));                      // e^{ik.R} = e^{2 pi i k_frac . n}
        }
    else
    {
        R.push_back(rvec3_t(0,0,0));
        ph.Append(dcmplx(1.0,0.0));
    }
    phase=ph.take();
    assert(!R.empty() && R.size()==phase.size());
}

// The local-PP sweeps' ABSOLUTE pair->level rule kappa (Ha per unit pair exponent; doc/GPWPlan.md 0e-PP
// step (a)): req = kappa*(alpha_i+alpha_j) bounds every pair's spectral tail by e^{-kappa/2} independent
// of the field's sharpness.  Default 30 Ha (e^{-15} -- CP2K's REL_CUTOFF 60 Ry); GPW_LOCALPP_RELCUTOFF
// overrides it for the self-convergence verification (kappa=60 must match to tolerance).  Read per call
// (not a static) so a test can A/B via setenv.
double LocalPPRelCutoff()
{
    double k=30.0;
    if (const char* e=std::getenv("GPW_LOCALPP_RELCUTOFF")) k=std::atof(e);
    return k;
}

// --- Real-space pseudopotential fields, replicated from the molecular PP_Local/PP_NonLocal terms (which live
//     Hamiltonian-side and so are out of reach from qcLattice_BS).  They are pure functions of the qcPseudo-
//     potential models + geometry, so the clean long-term home is qcPseudopotential (below both libraries);
//     the DRY-move is a deferred cleanup (see doc/GPWPlan.md).  Kept bit-identical to the term versions so a
//     Gaussian-in-a-box GPW PP matrix equals the finite molecular PP matrix. --------------------------------

// Real (tesseral) spherical harmonic Y_lm(rhat), unit-normalised on the sphere; the standard Cartesian forms
// (any orthonormal degree-l set is equivalent since only Sum_m |Y_lm><Y_lm| enters the KB projector).
double RealYlm(int l, int m, double x, double y, double z)
{
    using std::sqrt;
    const double pi=Pi;
    switch (l)
    {
    case 0: return 0.5/sqrt(pi);
    case 1:
        switch (m) { case -1: return sqrt(0.75/pi)*y; case 0: return sqrt(0.75/pi)*z; case 1: return sqrt(0.75/pi)*x; }
        break;
    case 2:
        switch (m)
        {
        case -2: return sqrt(15.0/(4*pi))*x*y;
        case -1: return sqrt(15.0/(4*pi))*y*z;
        case  0: return sqrt( 5.0/(16*pi))*(3*z*z-1.0);
        case  1: return sqrt(15.0/(4*pi))*x*z;
        case  2: return sqrt(15.0/(16*pi))*(x*x-y*y);
        }
        break;
    case 3:
        switch (m)
        {
        case -3: return sqrt( 35.0/(32*pi))*y*(3*x*x-y*y);
        case -2: return sqrt(105.0/(4 *pi))*x*y*z;
        case -1: return sqrt( 21.0/(32*pi))*y*(5*z*z-1.0);
        case  0: return sqrt(  7.0/(16*pi))*z*(5*z*z-3.0);
        case  1: return sqrt( 21.0/(32*pi))*x*(5*z*z-1.0);
        case  2: return sqrt(105.0/(16*pi))*z*(x*x-y*y);
        case  3: return sqrt( 35.0/(32*pi))*x*(x*x-3*y*y);
        }
        break;
    }
    assert(false && "RealYlm: unsupported (l,m)");
    return 0.0;
}

// One KB channel as a scalar field: beta_p(|r-R|) * Y_lm((r-R)^).  beta_p(0)=0 for l>=1 so the r=R angular
// singularity is harmless; guard the unit-vector division.
class BetaYlmField : public ScalarFunction<double>
{
    rvec3_t R; const Pseudopotential::SeparablePotential_R& v; int Z; size_t p; int l, m;
public:
    BetaYlmField(const rvec3_t& R_, const Pseudopotential::SeparablePotential_R& v_, int Z_, size_t p_, int l_, int m_)
        : R(R_), v(v_), Z(Z_), p(p_), l(l_), m(m_) {}
    double operator()(const rvec3_t& r) const override
    {
        rvec3_t d=r-R; double rr=norm(d);
        if (rr<1e-12) return l==0 ? v.BetaR(Z,p,0.0)*RealYlm(0,0,0,0,0) : 0.0;
        return v.BetaR(Z,p,rr) * RealYlm(l,m, d.x/rr, d.y/rr, d.z/rr);
    }
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}
};

// The radius beyond which the KB projector beta_p(r) is negligible (< 1e-10 of its peak).  beta_p is a GTH
// polynomial x Gaussian, compactly supported in practice; used to SCREEN the (image, mesh-point) projection
// loop in MakeSeparablePP.  Cheap (~1000 BetaR evals, once per projector).  1e-10 is far below any GPW anchor
// tolerance, so screening is numerically exact.
double BetaSupportRadius(const Pseudopotential::SeparablePotential_R& sep, int Z, size_t p)
{
    const double h=0.02, rmax=25.0, tol=1e-10;
    double peak=0.0;
    for (double r=0.0; r<=rmax; r+=h) peak=std::max(peak, std::fabs(sep.BetaR(Z,p,r)));
    double rsup=0.0;
    for (double r=0.0; r<=rmax; r+=h) if (std::fabs(sep.BetaR(Z,p,r))>tol*peak) rsup=r;
    return rsup+2.0*h;   // a small margin past the last significant radius
}

// --- Cartesian-monomial expansion of r^l Y_lm (the analytic-KB polynomial) --------------------------------
// The SOLID harmonic r^l Y_lm(rhat) is a homogeneous degree-l polynomial: Math::SphericalShell supplies the
// RAW monomial shape; the absolute scale is fixed NUMERICALLY against this file's own RealYlm table (evaluate
// both at a generic unit direction and take the ratio -- exact to machine precision since they are the same
// polynomial, and asserted at a second direction).  Guarantees the analytic KB uses the SAME Y_lm convention
// as the mesh path it replaces.
std::vector<qchem::Math::CartTerm> YlmCartesian(int l, int m)
{
    std::vector<qchem::Math::CartTerm> terms=qchem::Math::SphericalShell(l)[size_t(l+m)];
    auto poly=[&terms](double x, double y, double z)
    {
        double s=0.0;
        for (const auto& t : terms)
            s += t.c * std::pow(x,t.p.n) * std::pow(y,t.p.l) * std::pow(z,t.p.m);
        return s;
    };
    auto fix=[&](double x, double y, double z)                     // RealYlm(unit)*r^l / raw(x,y,z)
    {
        double r=std::sqrt(x*x+y*y+z*z);
        return RealYlm(l,m, x/r,y/r,z/r) * std::pow(r,l) / poly(x,y,z);
    };
    const double N=fix(0.3,0.5,0.7);                               // generic direction (no raw zeros, l<=3)
    assert(std::fabs(fix(0.9,-0.2,0.4)-N)<1e-10*std::fabs(N)+1e-14 &&
           "YlmCartesian: raw solid harmonic is not proportional to RealYlm (convention drift)");
    for (auto& t : terms) t.c*=N;
    return terms;
}

// Multiply a Cartesian polynomial by (x^2+y^2+z^2)^n (the r^{2n} factor of a higher radial projector),
// combining duplicate monomials.
std::vector<qchem::Math::CartTerm> MultiplyR2(std::vector<qchem::Math::CartTerm> terms, int n)
{
    using qchem::Math::Monomial;
    for (int k=0;k<n;k++)
    {
        std::map<Monomial,double> acc;
        for (const auto& t : terms)
            for (int ax=0; ax<3; ax++)
            {
                Monomial p=t.p; p[ax]+=2;
                acc[p]+=t.c;
            }
        terms.clear();
        for (const auto& [p,c] : acc) terms.push_back({p,c});
    }
    return terms;
}
} //anon

GPW_Evaluator::GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                             double densityEcut, const rvec3_t& kFrac, bool kIsReal, bool homeCellOnly,
                             double cutoffFactor, RasterPolicy raster, double ladderFactor,
                             RasterFields rasterFields)
    : itsMol(std::move(mol))
    , itsHomeOnly(homeCellOnly)
    , itsk(kFrac)
    , itsTRIM(kIsReal)
    , itsCell(cell)
    , itsCutoffFactor(cutoffFactor)
    , itsLadderFactor(ladderFactor)
    , itsRaster(raster)
{
    // The single orbital block of the (raw, no-SALC) molecular Gaussian basis.
    const BasisSet::Real_OIBS* only=nullptr;
    for (auto ibs : itsMol->Iterate<BasisSet::Real_OIBS>()) { assert(!only); only=ibs; }
    if (!only) throw std::runtime_error("GPW_Evaluator: no orbital block in the molecular basis");
    itsOrb=only;
    itsN  =itsOrb->GetNumFunctions();

    // Its periodic-1E capability, reached by an abstract->abstract cross-cast (a molecular Gaussian basis
    // realises Molecule::LatticeSum1E; anything else is a usage error).
    itsLat=dynamic_cast<const Molecule::LatticeSum1E*>(itsOrb);
    if (!itsLat) throw std::runtime_error(
        "GPW_Evaluator: the orbital basis is not a molecular Gaussian basis (no Molecule::LatticeSum1E)");

    // THERE IS NO CUT (user pin, doc/GPWPlan.md): every lattice sum is an eps-CONVERGED SERIES enumerated
    // INSIDE the molecular seam per shell pair (1E matrices, analytic KB, collocation) -- no radius
    // parameter exists.  The ONE remaining explicit image list is the INTERNAL Bloch-orbital set for
    // Eval/EvalGradient + the mesh-path KB quadrature, DERIVED from the same eps-screen: a single orbital
    // reaches sqrt(-ln eps/alpha_min) and its centre sits within a cell span of any evaluation point, so
    // reach + span covers every image the screen keeps (screening then prunes it sparse per point).  The
    // home-only MODE (the finite-molecule configuration) keeps just the origin.
    // VERIFICATION INSTRUMENT (doc/GPWPlan 0.5(a); the GPW_MGRID_ECUTS precedent): GPW_RASTER_POLICY=ball
    // flips EVERY grid this block builds (density grid + ladder levels) to the BallOnly raster -- the A/B
    // knob for the raster-policy calibration.  The shipped surface stays AliasFree until the A/B verdict
    // promotes the enum onto the factory signatures (or closes the item as policy-justified).
    if (const char* rp=std::getenv("GPW_RASTER_POLICY"))
        itsRaster = std::string(rp)=="ball" ? RasterPolicy::BallOnly : RasterPolicy::AliasFree;

    // HartreeOnly routing beta (needs itsLat's alpha_max, so set here, not in the init list).  MEASURED
    // (NaF SR2 Becke, 2026-07-31): beta=0 (pure pair-bandwidth) DIVERGES (+904 Ha, low-G slosh) -- an
    // under-resolved diffuse pair's FFT folds high-G back into the kept ball, corrupting rho-tilde's low-G
    // (the charge-transfer mode), so the 2/3*alpha_max floor was protecting the DENSITY too, not only
    // V_xc.  The HartreeOnly floor is therefore a FRACTION of alpha_max, sweepable via GPW_RELFIELDSHARP
    // (the calibration instrument; pin by convergence + rho_lost/N, not wall-clock).
    if (rasterFields==RasterFields::HartreeOnly)
    {
        const char* s=std::getenv("GPW_RELFIELDSHARP");
        const double frac = s ? std::atof(s) : 1.0/3.0;
        itsRelFieldSharp = frac*itsLat->MaxExponent();
    }
    itsMaxReach=std::sqrt(-std::log(1e-10)/itsLat->MinExponent());
    itsCellCtr =cell.ToCartesian(rvec3_t(0.5,0.5,0.5));
    itsCellRad =0.0;
    for (int cx=0;cx<=1;cx++) for (int cy=0;cy<=1;cy++) for (int cz=0;cz<=1;cz++)
    {
        rvec3_t corner=cell.ToCartesian(rvec3_t(double(cx),double(cy),double(cz)));
        itsCellRad=std::max(itsCellRad, norm(corner-itsCellCtr));
    }
    // Image-list radius: Eval's keep test is |r - Rc - ctr| <= cellRad+maxReach with r in (or near) the
    // cell, so the exact bound is |Rc| <= maxReach + 2*cellRad.  For an OBLIQUE cell (MnO rhombohedral)
    // 2*cellRad (the cell DIAMETER) exceeds 2*maxCellEdge, so the historical formula under-enumerated;
    // keep the max of both so every previously-enumerated case is unchanged.
    BuildImages(cell, itsHomeOnly ? 0.0 : std::max(2.0*itsMaxReach+2.0*cell.GetMaximumCellEdge(),
                                                   itsMaxReach+2.0*itsCellRad), itsk, itsTRIM, itsRc, itsPhaseC);

    // The DFT tier's density/collocation grid: GPW's ONLY grid cutoff (no orbital/wavefunction cutoff -- the
    // Gaussians are analytic).  IMPORTANT: densityEcut is a DENSITY-scale quantity: the sharpest feature is
    // the PRODUCT of the two tightest primitives (a Gaussian of exponent 2*alpha_max), and cutoffFactor=2
    // (the default) resolves the density AT its own exponent.  The pre-(f2) default of 8 was compensation
    // for the ball-projected XC feed's Gibbs lobes (see the ctor doc + doc/GPWPlan 0.5(f)); with the raw
    // collocated feed the measured NaF curve is a sub-2-mHa plateau from C~1.5 up.
    //   densityEcut < 0 : AUTOMATIC (recommended) -- floor = cutoffFactor*alpha_max, from the basis (no user Ha).
    //   densityEcut = 0 : DFT tier OFF (1E-only; no grid).
    //   densityEcut > 0 : EXPLICIT -- honoured as given, but WARN on cerr if below the floor (the caller insisted
    //                     on an under-resolved grid: charge leaks off-grid -- we don't hide it, but we don't
    //                     silently override the explicit choice either).  For SIPP alpha_max=2 the floor is 4,
    //                     below every committed Si anchor's densityEcut (>=6), so those are unchanged.
    if (densityEcut!=0.0)
    {
        const double aMax  = itsLat->MaxExponent();
        const double floor = cutoffFactor * aMax;          // cutoffFactor*alpha_max (DENSITY scale)
        double ecut = densityEcut;
        if (densityEcut < 0.0)                              // AUTOMATIC: pick the basis-derived floor
            ecut = floor;
        else if (densityEcut < floor)                      // EXPLICIT but under-resolved: warn, honour as given
            std::cerr << "[GPW] WARNING: densityEcut=" << densityEcut << " < " << cutoffFactor
                      << "*alpha_max=" << floor << " (basis alpha_max=" << aMax << "): the density grid "
                      << "under-resolves the basis -- charge will leak off-grid (integral rho_grid < N). Prefer "
                      << "densityEcut>=" << floor << ", or the automatic default (densityEcut<0)." << std::endl;
        itsFFT_R_G_Grids=std::make_shared<const PW_Grid_Evaluator>(
                    ReciprocalLattice(cell.MakeReciprocalCell()), rvec3_t(0,0,0), ecut, itsRaster);
    }
}

// Out of line (the class is polymorphic and the dtor is the key function).  It used to hand this block's
// ladder-shaped collocation STREAMS back to a global point budget (doc/GPWPlan.md 0.5(b)); the value cache
// is gone (doc/CollocationRewritePlan.md step 7) and its replacement -- the (shell pair, offset) task list --
// is a few hundred kB of pure geometry with nothing to refund, so there is nothing left to release.
GPW_Evaluator::~GPW_Evaluator() = default;

// The Bloch phase closure the analytic kernels call back per screened cross-cell offset: e^{2 pi i k_frac.n}
// (n the INTEGER cell index -- convention-safe, same convention as BuildImages).  Gamma: the constant 1.
Molecule::LatticeSum1E::cellphase_t GPW_Evaluator::CellPhase() const
{
    const rvec3_t k=itsk;
    const bool trim=itsTRIM;                         // the irrep's exact fact; phases exactly +/-1 there (BlochPhase)
    return [k,trim](const ivec3_t& n)->dcmplx { return BlochPhase(k,n,trim); };
}

// The matrix-free density->rho-tilde / ->V_H map (the Projector3<dcmplx>::apply realization -- the CP2K analytic
// collocation, doc/GPWPlan.md S0).  The density is collocated ANALYTICALLY per (pair, screened cross-cell
// offset) on compact exp-tail boxes, modulo-wrapped, each pair on the coarsest REL_CUTOFF level resolving its
// product exponent (LatticeSum1E::CollocateDensity); each level is FFT'd and rho-tilde is combined NESTED in
// G-space (a coarse level contributes over its own {G} -- beyond it the pair's spectrum is below tolerance by
// construction).  rho-tilde(G)=FFT[rho]/Npts (the grid quadrature's 1/Omega and Omega/Npts cancel).  For
// Coulomb the diagonal Poisson kernel 4pi/|G|^2 is folded in (V_H = 4pi rho-tilde/G^2; G=0 -> 0).  The closure
// keeps the level grids + molecular basis alive (captured shared_ptrs) since it lives in the framework-cached
// Projector3<dcmplx>.
std::function<ΔG_Map(const chmat_t&)> GPW_Evaluator::MakeCollocator(bool coulomb, std::shared_ptr<const PW_Grid_Evaluator> grid) const
{
    // The ladder derives from the REQUESTED fit grid (not the block's own itsFFT_R_G_Grids); the closure captures it.
    std::vector<std::shared_ptr<const PW_Grid_Evaluator>> levels;
    std::vector<ivec3_t> N_L;
    std::vector<double>  ecut_L;
    size_t nBase=0;
    BuildLevels(grid, levels, N_L, ecut_L, nBase);
    const UnitCell A = itsCell;                             // the direct cell (stored; see the member note)
    auto mol    = itsMol;                                   // shared_ptr -> keep the molecular basis (lat) alive
    const Molecule::LatticeSum1E* lat = itsLat;
    auto phase  = CellPhase();
    ReciprocalLattice recip = grid->Recip();
    if (!itsCollocMemo) itsCollocMemo=std::make_shared<CollocMemo>();
    auto memo   = itsCollocMemo;                            // ONE memo across the Coulomb + overlap closures
    return [A,levels,N_L,ecut_L,mol,lat,phase,recip,coulomb,memo,relFS=itsRelFieldSharp](const chmat_t& D) -> ΔG_Map
    {
        // ⚠ CALL CENSUS (2026-08-28).  These four closures are the ONLY consumers of the analytic
        // collocate/integrate pair, so bucketing their bodies gives the ledger a per-SITE call COUNT for
        // free -- and, because the inner CollocateDensity/IntegratePotential buckets are exclusive
        // CHILDREN, the site bucket's own time is exactly the FFT + spectral-transfer work that was
        // previously unattributed.  Two answers from one label; see doc/OpenWork.md's call-census note.
        qchem::report::Timed timer("scf: rho ball closure (FFT + G-combine)");
        // Replay a density this evaluator has ALREADY collocated on this ladder (the sibling tensor's call,
        // or the other spin channel's).  Exact equality or recompute -- see CollocMemo::Lookup.
        std::vector<rvec_t> rho;
        if (!memo->Lookup(D, ecut_L, rho))
        {
            rho = lat->CollocateDensity(D, phase, A, N_L, ecut_L, relFS);
            memo->Store(D, rho, ecut_L);
        }
        else { qchem::report::Timed t("scf: rho memo replay (D seen before)"); }
        ΔG_Map out;
        for (size_t l=0; l<levels.size(); l++)
        {
            cvec_t rhoTilde = levels[l]->ForwardFFT(rho[l]);            // this level's rho-tilde(G)
            for (const ivec3_t& dm : levels[l]->Gs())                   // nested combine over the level's own {G}
                out[dm] += levels[l]->GridCoeff(rhoTilde, dm);
        }
        if (coulomb)
            for (auto& [dm,c] : out) c *= recip.CoulombKernel(dm);
        return out;
    };
}
// Coulomb 3-centre tensor: the matrix-free analytic collocation map with the Poisson kernel folded in.  The
// columns still list the fit-basis {G} (the fine grid's); the per-column kernel is inside `apply` (which
// short-circuits Contract), so `kernel` stays empty.
// Over the REQUESTED fit grid: columns are grid's {G}, the collocation ladder derives from grid.  The no-arg
// overloads below are the block's-own-grid convenience (Repulsion3CTensor(itsFFT_R_G_Grids)).
Projector3<dcmplx> GPW_Evaluator::Repulsion3CTensor(std::shared_ptr<const PW_Grid_Evaluator> grid) const
{
    Projector3<dcmplx> g;
    g.volume=grid->Volume();
    for (const ivec3_t& dm : grid->Gs()) g.columns.push_back({dm,{}});
    g.apply       =MakeCollocator(/*coulomb*/true, grid);   // forward: D -> V_H (Coulomb baked)
    g.applyAdjoint=MakeIntegrator(grid);                    // backward: field -> <i|f|j> (overlap; kernel fwd-only)
    return g;
}
// Overlap 3-centre tensor: the same analytic map, no kernel (the density's own rho-tilde).
Projector3<dcmplx> GPW_Evaluator::Overlap3CTensor(std::shared_ptr<const PW_Grid_Evaluator> grid) const
{
    Projector3<dcmplx> g;
    g.volume=grid->Volume();
    for (const ivec3_t& dm : grid->Gs()) g.columns.push_back({dm,{}});
    g.apply       =MakeCollocator(/*coulomb*/false, grid);  // forward: D -> rho-tilde
    g.applyAdjoint=MakeIntegrator(grid);                    // backward: v_xc-tilde -> <i|v_xc|j>
    g.applyRaw       =MakeRawCollocator(grid);              // 0.5(f2): D -> rho_DM(r) raw (XC feed)
    g.applyRawAdjoint=MakeRawIntegrator(grid);              //          v(r) -> <i|v|j> (exact transpose)
    return g;
}

// The Projector3<dcmplx> BACKWARD realization (the adjoint of MakeCollocator): a SELF-CONTAINED integrate-back closure over
// the fit \a grid's REL_CUTOFF ladder -- <i|f|j> = analytic per-pair gather of the field f (LatticeSum1E::
// IntegratePotential), the exact adjoint of the collocation on the SAME grid.  So Overlap3C/Repulsion3C(fit)
// carry BOTH directions on the fit grid; the KS matrix is built by ContractAdjoint (doc/GPWPlan §0e step2).
// (OverlapMatrix below is the block's-own-itsFFT_R_G_Grids sibling that MakeLocalPP + the legacy MakeOverlap still use.)
std::function<chmat_t(const std::function<dcmplx(const ivec3_t&)>&)>
GPW_Evaluator::MakeIntegrator(std::shared_ptr<const PW_Grid_Evaluator> grid) const
{
    std::vector<std::shared_ptr<const PW_Grid_Evaluator>> levels;
    std::vector<ivec3_t> N_L;
    std::vector<double>  ecut_L;
    size_t nBase=0;
    BuildLevels(grid, levels, N_L, ecut_L, nBase);          // the fit grid's ladder
    const UnitCell A = itsCell;
    auto mol   = itsMol;                                    // keep the molecular basis (lat) alive in the closure
    const Molecule::LatticeSum1E* lat = itsLat;
    auto phase = CellPhase();
    if (!itsCollocMemo) itsCollocMemo=std::make_shared<CollocMemo>();
    if (!itsGatherMemo) itsGatherMemo=std::make_shared<GatherMemo>();
    auto memo  = itsCollocMemo;
    auto gmemo = itsGatherMemo;                             // shared D-aware screen (adjoint of the collocation)
    return [A,levels,N_L,ecut_L,mol,lat,phase,memo,gmemo,relFS=itsRelFieldSharp](const std::function<dcmplx(const ivec3_t&)>& Vtilde) -> chmat_t
    {
        qchem::report::Timed timer("scf: h ball closure (level restrict + inverse FFT)");
        const size_t K=levels.size();
        std::vector<rvec_t> V_L(K);
        for (size_t L=0;L<K;L++)
        {
            ΔG_Map vmapL;
            for (const ivec3_t& dm : levels[L]->Gs()) vmapL[dm]=Vtilde(dm);   // restrict to level L's {G}
            V_L[L]=levels[L]->RhoOnGrid(vmapL);
        }
        const chmat_t* screenD = (memo && memo->valid) ? &memo->Dscr : nullptr;   // the UNION screen (all channels; see CollocMemo)
        chmat_t h;
        if (gmemo->Lookup(V_L, screenD, h)) { qchem::report::Timed t("scf: h memo replay (V+screen seen)"); return h; }
        h = lat->IntegratePotential(V_L, phase, A, N_L, ecut_L, 0.0, screenD, 0.0, relFS);   // 0 = the relative rule (adjoint-paired with collocation)
        gmemo->Store(V_L, screenD, h);
        return h;
    };
}
// doc/GPWPlan 0.5(f2): spectral transfer between two FFT rasters, in ForwardFFT coefficient layout.
// Copies every mode in the ALIAS-FREE band of BOTH boxes (2|m| < N per axis on each side -- the same strict
// band RhoOnGrid/GridCoeff speak) and drops the rest: source->finer = exact zero-pad upsampling of the
// trigonometric interpolant; source->coarser = band truncation.  The two directions are TRANSPOSES of each
// other, which is what makes the raw-XC H_xc the exact derivative of the raw-XC energy.
static void TransferBand(const cvec_t& src, const ivec3_t& Ns, cvec_t& dst, const ivec3_t& Nd)
{
    const int mx=(Ns.x-1)/2, my=(Ns.y-1)/2, mz=(Ns.z-1)/2;            // strict alias-free band of the source
    const int dx=(Nd.x-1)/2, dy=(Nd.y-1)/2, dz=(Nd.z-1)/2;            //   and of the destination
    auto w=[](int m, int N) { return size_t(((m%N)+N)%N); };
    for (int ix=-std::min(mx,dx); ix<=std::min(mx,dx); ix++)
        for (int iy=-std::min(my,dy); iy<=std::min(my,dy); iy++)
            for (int iz=-std::min(mz,dz); iz<=std::min(mz,dz); iz++)
                dst[(w(ix,Nd.x)*Nd.y+w(iy,Nd.y))*Nd.z+w(iz,Nd.z)] +=
                src[(w(ix,Ns.x)*Ns.y+w(iy,Ns.y))*Ns.z+w(iz,Ns.z)];
}

// --- CollocMemo: the density replay cache (doc/OpenWork.md, the 2026-08-28 call census) ---------------
// ONE copy of the match/store logic.  It used to be written out twice -- once in MakeCollocator, once in
// MakeRawCollocator with the comment "the MakeCollocator memo logic, verbatim" -- and this tree has already
// paid once for a second copy of a screening rule that drifted (the integrate-back Re[D conj(phase)] defect,
// doc/CollocationRewritePlan.md step 7).  Two call sites, one definition.
size_t GPW_Evaluator::CollocMemo::Depth()
{
    static const size_t d=[]{ const char* s=std::getenv("GPW_COLLOC_MEMO"); return s ? size_t(std::atoi(s)) : 4u; }();
    return d;
}
namespace {
//! EXACT equality on the upper triangle -- never a relaxed compare.  A replay must be bit-identical to the
//! collocation it stands in for, or the memo becomes a silent accuracy knob.
bool SameDensity(const chmat_t& a, const chmat_t& b)
{
    if (a.rows()!=b.rows()) return false;
    for (size_t i=0;i<a.rows();i++)
        for (size_t j=i;j<a.columns();j++) if (a(i,j)!=b(i,j)) return false;
    return true;
}
} // anon
bool GPW_Evaluator::CollocMemo::Lookup(const chmat_t& d, const std::vector<double>& ec,
                                       std::vector<rvec_t>& rhoOut)
{
    // The ladder is part of the key: the SAME D on a DIFFERENT ladder has a different rho (the cross-block
    // seed hands one density to two grids), and replaying the wrong grid is the one way this can be wrong.
    if (valid && ecut==ec && SameDensity(D,d)) { rhoOut=rho; return true; }
    for (size_t k=0;k<past.size();k++)
        if (past[k].ecut==ec && SameDensity(past[k].D,d))
        {
            Entry hit=std::move(past[k]);                  // PROMOTE: `D` must stay "the most recent"
            past.erase(past.begin()+long(k));
            if (valid) past.insert(past.begin(), Entry{D, rho, ecut});
            D=hit.D; rho=hit.rho; ecut=hit.ecut; valid=true;
            if (past.size()>Depth()) past.resize(Depth());
            rhoOut=rho;
            return true;
        }
    return false;
}
void GPW_Evaluator::CollocMemo::Store(const chmat_t& d, const std::vector<rvec_t>& r,
                                      const std::vector<double>& ec)
{
    // Dscr = the UNION screen: reset on a new ladder/shape, else WIDEN to cover this channel's active set
    // too (a polarized run collocates D_up then D_dn through the same memo).  Widen-only is a PIN -- see
    // the CollocMemo interface note; an epoch reset was tried and rejected.
    const bool reset = !valid || Dscr.rows()!=d.rows() || ecut!=ec;
    if (reset) Dscr=chmat_t(d.rows());
    for (size_t i=0;i<d.rows();i++)
        for (size_t j=i;j<d.columns();j++)
            Dscr(i,j) = reset ? std::abs(dcmplx(d(i,j)))
                              : std::max(std::real(dcmplx(Dscr(i,j))), std::abs(dcmplx(d(i,j))));
    if (valid) past.insert(past.begin(), Entry{D, rho, ecut});   // retire the previous most-recent
    if (past.size()>Depth()) past.resize(Depth());
    D=d; rho=r; ecut=ec; valid=true;
}

bool GPW_Evaluator::GatherMemo::Lookup(const std::vector<rvec_t>& V_L, const chmat_t* screen,
                                       chmat_t& hOut) const
{
    auto sameField=[&](const std::vector<rvec_t>& a)->bool
    {
        if (a.size()!=V_L.size()) return false;
        for (size_t l=0;l<a.size();l++)
        {
            if (a[l].size()!=V_L[l].size()) return false;
            for (size_t k=0,m=a[l].size(); k<m; k++) if (a[l][k]!=V_L[l][k]) return false;   // exact
        }
        return true;
    };
    auto sameScreen=[&](const Entry& e)->bool
    {
        if (e.hasScreen != (screen!=nullptr)) return false;
        if (!screen) return true;
        if (e.screen.rows()!=screen->rows()) return false;
        for (size_t i=0;i<screen->rows();i++)
            for (size_t j=i;j<screen->columns();j++) if (e.screen(i,j)!=(*screen)(i,j)) return false;
        return true;
    };
    for (const Entry& e : entries)
        if (sameScreen(e) && sameField(e.V_L)) { hOut=e.h; return true; }
    return false;
}
void GPW_Evaluator::GatherMemo::Store(const std::vector<rvec_t>& V_L, const chmat_t* screen, const chmat_t& h)
{
    Entry e;
    e.V_L=V_L; e.h=h; e.hasScreen=(screen!=nullptr);
    if (screen) e.screen=*screen;
    entries.insert(entries.begin(), std::move(e));
    if (entries.size()>CollocMemo::Depth()) entries.resize(CollocMemo::Depth());
}

// D -> rho_DM(r) on grid's raster: level 0 (== grid) RAW, every other level spectrally transferred in.
std::function<rvec_t(const chmat_t&)> GPW_Evaluator::MakeRawCollocator(std::shared_ptr<const PW_Grid_Evaluator> grid) const
{
    std::vector<std::shared_ptr<const PW_Grid_Evaluator>> levels;
    std::vector<ivec3_t> N_L;
    std::vector<double>  ecut_L;
    size_t nBase=0;
    BuildLevels(grid, levels, N_L, ecut_L, nBase);
    const UnitCell A = itsCell;
    auto mol   = itsMol;
    const Molecule::LatticeSum1E* lat = itsLat;
    auto phase = CellPhase();
    if (!itsCollocMemo) itsCollocMemo=std::make_shared<CollocMemo>();
    if (!itsGatherMemo) itsGatherMemo=std::make_shared<GatherMemo>();
    auto memo  = itsCollocMemo;
    auto gmemo = itsGatherMemo;
    return [A,levels,N_L,ecut_L,mol,lat,phase,memo,relFS=itsRelFieldSharp](const chmat_t& D) -> rvec_t
    {
        qchem::report::Timed timer("scf: rho raw closure (spectral transfer)");
        std::vector<rvec_t> rho;                            // ONE memo, ONE definition (see CollocMemo::Lookup)
        if (!memo->Lookup(D, ecut_L, rho))
        {
            rho = lat->CollocateDensity(D, phase, A, N_L, ecut_L, relFS);
            memo->Store(D, rho, ecut_L);
        }
        else { qchem::report::Timed t("scf: rho memo replay (D seen before)"); }
        const ivec3_t NT=N_L[0];                            // the integration raster (level 0 == grid)
        cvec_t acc(size_t(NT.x)*NT.y*NT.z, dcmplx(0.0));
        for (size_t l=1; l<levels.size(); l++)              // coarse levels up, top rung band-dropped down
            TransferBand(levels[l]->ForwardFFT(rho[l]), N_L[l], acc, NT);
        rvec_t out=levels[0]->BackwardFFT(acc);
        out+=rho[0];                                        // the finest level RAW: no FFT round trip at all
        return out;
    };
}

// The exact transpose: v(r) on grid's raster -> <i|v|j>.  V_0 = v identity (adjoint of the raw keep);
// V_L = BackwardFFT_L(TransferBand(ForwardFFT_T(v))) -- the ForwardFFT 1/Npts against IntegratePotential's
// internal Omega/Npts(L) weight leaves NO stray factors (derivation in the interface doc).
std::function<chmat_t(const rvec_t&)> GPW_Evaluator::MakeRawIntegrator(std::shared_ptr<const PW_Grid_Evaluator> grid) const
{
    std::vector<std::shared_ptr<const PW_Grid_Evaluator>> levels;
    std::vector<ivec3_t> N_L;
    std::vector<double>  ecut_L;
    size_t nBase=0;
    BuildLevels(grid, levels, N_L, ecut_L, nBase);
    const UnitCell A = itsCell;
    auto mol   = itsMol;
    const Molecule::LatticeSum1E* lat = itsLat;
    auto phase = CellPhase();
    if (!itsCollocMemo) itsCollocMemo=std::make_shared<CollocMemo>();
    if (!itsGatherMemo) itsGatherMemo=std::make_shared<GatherMemo>();
    auto memo  = itsCollocMemo;
    auto gmemo = itsGatherMemo;
    return [A,levels,N_L,ecut_L,mol,lat,phase,memo,gmemo,relFS=itsRelFieldSharp](const rvec_t& v) -> chmat_t
    {
        qchem::report::Timed timer("scf: h raw closure (spectral transfer)");
        const size_t K=levels.size();
        const ivec3_t NT=N_L[0];
        assert(v.size()==size_t(NT.x)*NT.y*NT.z);
        cvec_t vt=levels[0]->ForwardFFT(v);
        std::vector<rvec_t> V_L(K);
        V_L[0]=v;
        for (size_t l=1; l<K; l++)
        {
            cvec_t ctL(size_t(N_L[l].x)*N_L[l].y*N_L[l].z, dcmplx(0.0));
            TransferBand(vt, NT, ctL, N_L[l]);
            V_L[l]=levels[l]->BackwardFFT(ctL);
        }
        const chmat_t* screenD = (memo && memo->valid) ? &memo->Dscr : nullptr;   // the UNION screen (all channels; see CollocMemo)
        chmat_t h;
        if (gmemo->Lookup(V_L, screenD, h)) { qchem::report::Timed t("scf: h memo replay (V+screen seen)"); return h; }
        h = lat->IntegratePotential(V_L, phase, A, N_L, ecut_L, 0.0, screenD, 0.0, relFS);
        gmemo->Store(V_L, screenD, h);
        return h;
    };
}

Projector3<dcmplx> GPW_Evaluator::Repulsion3CTensor() const {return Repulsion3CTensor(itsFFT_R_G_Grids);}
Projector3<dcmplx> GPW_Evaluator::Overlap3CTensor  () const {return Overlap3CTensor  (itsFFT_R_G_Grids);}

// The potential->KS bridge (the applyAdjoint realization / isPW_DFT_Evaluator) -- the ANALYTIC integrate-back,
// the exact adjoint of MakeCollocator's density collocation.  Restrict Vtilde to each REL_CUTOFF level's own
// {G} (a SPECTRAL low-pass -- no ringing), inverse-FFT to that level's grid, and let the molecular side gather
// each orbital pair analytically on ITS level (LatticeSum1E::IntegratePotential -- same boxes, same wrap, same
// level assignment as the density side, so <collocate(D),V> == <D,integrate(V)> to machine precision).  Only V
// is ever sampled, so even the SHARP local PP integrates accurately against a coarse-level diffuse pair (the
// pair's own bandwidth bounds what it senses of V) -- MakeLocalPP routes through here too.
chmat_t GPW_Evaluator::OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    EnsureLevels();
    const UnitCell& A=itsCell;                              // the direct cell (stored; see the member note)
    const size_t K=itsLevels.size();
    std::vector<rvec_t> V_L(K);
    for (size_t L=0;L<K;L++)
    {
        ΔG_Map vmapL;
        for (const ivec3_t& dm : itsLevels[L]->Gs()) vmapL[dm]=Vtilde(dm);   // restrict to level L's {G} (low-pass)
        V_L[L]=itsLevels[L]->RhoOnGrid(vmapL);
    }
    // D-AWARE integrate-back: the KS fields integrated here derive from the density last collocated by this
    // evaluator (the CollocMemo holds the iteration's D), so pass it as the seam's density screen -- both
    // directions then keep the IDENTICAL active set (adjoint exact on the shared truncated operator) and the
    // sweep skips every term the density cannot resolve.  Before the first collocation (or for a field
    // unrelated to a density, e.g. the unit-field gates) the memo is empty -> complete sweep.
    const chmat_t* screenD = (itsCollocMemo && itsCollocMemo->valid) ? &itsCollocMemo->D : nullptr;
    return itsLat->IntegratePotential(V_L, CellPhase(), A, itsLevelN, itsLevelEcut, 0.0, screenD, 0.0, itsRelFieldSharp);   // 0 = the relative rule (adjoint-paired with collocation)
}

// The REL_CUTOFF multi-grid density-grid ladder: the fine grid (L=0, reused) plus coarser grids each a factor
// 4 lower in Ecut, down to the level that still resolves the most-diffuse orbital product (exponent 2*alpha_min).
// Built once (geometry-fixed).  Because the collocation/integrate-back are ANALYTIC (only V is ever sampled), a
// diffuse pair on its matched coarse level is accurate -- the sampling multigrid's aliasing (and its depth cap)
// are gone.
// The block's OWN ladder (from itsFFT_R_G_Grids), cached for OverlapMatrix (the integrate-back) + MakeLocalPP.  The
// collocator instead builds its ladder from the REQUESTED fit grid (MakeCollocator -> BuildLevels(grid)).
void GPW_Evaluator::EnsureLevels() const
{
    if (!itsLevels.empty()) return;
    assert(itsFFT_R_G_Grids && "GPW_Evaluator: the DFT tier requires the density grid (densityEcut!=0)");
    BuildLevels(itsFFT_R_G_Grids, itsLevels, itsLevelN, itsLevelEcut, itsNBaseLevels);
}

// The grid-parameterized ladder core: level [0]==\a grid (the fine grid), coarser levels a factor 4 down to the
// most-diffuse pair floor, plus the top completion rung -- built from \a grid, NOT the block's itsFFT_R_G_Grids, so the
// collocator honours the requested fit grid.  (Orbital-exponent + cell state -- MaxExponent/MinExponent/itsCell/
// itsCutoffFactor/RelCutoffSafety -- stays the block's; only the GRID is parameterized.)
void GPW_Evaluator::BuildLevels(std::shared_ptr<const PW_Grid_Evaluator> grid,
                                std::vector<std::shared_ptr<const PW_Grid_Evaluator>>& levels,
                                std::vector<ivec3_t>& levelN, std::vector<double>& levelEcut, size_t& nBaseLevels) const
{
    assert(grid && "GPW_Evaluator::BuildLevels requires a density grid (densityEcut!=0)");
    // GRID-MATCHING OVERRIDE (doc/GPWPlan §0e "grid-matched CP2K validation"; verification instrument, not an
    // interface): GPW_MGRID_ECUTS="53.33,17.78,5.926" (Ha, comma-separated, descending) replaces the factor-4
    // sub-level progression AND the top completion rung with the EXPLICIT list -- level 0 stays \a grid (the
    // reference), the listed cutoffs follow.  CP2K's ladder is CUTOFF/3^(i-1) (progression_factor default 3),
    // which the factor-4 default cannot reproduce.  Pair->level assignment matching is the separate
    // GPW_RELCUTOFF knob (molecular side).
    if (const char* s=std::getenv("GPW_MGRID_ECUTS"))
    {
        levels.push_back(grid);
        for (std::string list(s); !list.empty();)
        {
            const size_t c=list.find(',');
            const double e=std::atof(list.substr(0,c).c_str());
            list = c==std::string::npos ? std::string() : list.substr(c+1);
            assert(e>0.0 && "GPW_MGRID_ECUTS: sub-level cutoffs must be positive");
            // A listed cutoff at/above the reference is SKIPPED (warn): the knob is process-wide, so a
            // coarser block in the same run (e.g. the grid-continuation SEED stage) keeps a valid ladder.
            if (e>=grid->Ecut())
            {
                std::cerr<<"[GPW] GPW_MGRID_ECUTS: skipping sub-level "<<e<<" >= reference "<<grid->Ecut()<<std::endl;
                continue;
            }
            levels.push_back(std::make_shared<const PW_Grid_Evaluator>(grid->Recip(), rvec3_t(0,0,0), e, itsRaster));
        }
        nBaseLevels=levels.size();          // no top rung: the explicit list IS the whole ladder
        for (const auto& g : levels)
        {
            levelN.push_back(g->FFTGrid());
            levelEcut.push_back(g->Ecut());
        }
        return;
    }
    const double efine=grid->Ecut();
    const double amax=itsLat->MaxExponent(), amin=itsLat->MinExponent();
    const double ecoarse=efine*amin/amax;             // resolves the most-diffuse pair product (exponent 2*amin)
    levels.push_back(grid);                           // L=0: the fine grid, reused
    // AUTOMATIC LADDER DEPTH (no user knob -- doc/GPWPlan §1 "diffuse at will").  Factor-4 steps from efine down
    // to the diffuse floor ecoarse span the alpha_max/alpha_min ratio, so the count is ~log4(amax/amin) coarse
    // levels -- the "e-folding" heuristic: Si's narrow exponent range self-selects ~1-2, ionic NaF's 1333x range
    // self-selects several.  Each level's own Ecut resolves the diffuse band the REL_CUTOFF rule assigns it.
    const double f=itsLadderFactor>1.0 ? itsLadderFactor : 4.0;   // progression factor (GPWParams.ladderFactor)
    double e=efine;
    while (e/f>=ecoarse)                              // factor-f coarsening down to the diffuse floor
    {
        e/=f;
        auto g=std::make_shared<const PW_Grid_Evaluator>(grid->Recip(), rvec3_t(0,0,0), e, itsRaster);
        // DEGENERACY FLOOR (min grid points): the ONLY failure mode of a coarse level is a cell too small to
        // hold a coarser grid -- an FFT grid saturating at N~1-2 no longer scales like sqrt(Ecut) and loses
        // charge (a sigma~2 pair on h=7 lost percent-level charge -- the failed crystal gate).  Keep the
        // properly-coarse N>=kMinLevelN grids the ultra-diffuse pairs want; stop below (those pairs then use the
        // coarsest surviving level, finer than they need).  Replaces an amax-driven h^2*pmax check whose value
        // ~= pi^2*amax/efine was ~LEVEL-INDEPENDENT, so it broke on the FIRST coarse level for every real basis
        // -- collapsing the ladder to L0+rung (the "we used to have 4 levels" regression).
        static constexpr int kMinLevelN=3;
        const ivec3_t Ng=g->FFTGrid();
        if (std::min(std::min(Ng.x,Ng.y),Ng.z) < kMinLevelN) break;
        levels.push_back(g);
    }
    nBaseLevels=levels.size();
    // TOP COMPLETION RUNG (doc/GPWPlan.md 0b').  The pair->level rule demands ecut >= RelCutoffSafety()
    // * efine * (a_i+a_j)/(2 a_max), whose maximum (the a_max+a_max pair) is RelCutoffSafety()*efine --
    // ABOVE the reference grid, so without this rung the sharpest pairs sit on the reference carrying an
    // e^{-ecut/2p} energy-tail (an ENERGY effect ONLY: collocated CHARGE is ball-truncation-immune, gated
    // by GPW.SharpestPairChargeConservation).  One rung at exactly that cutoff makes every smooth-path
    // assignment satisfiable.  Appended LAST so ecut_L[0] stays the resolution reference (the requirement
    // must not stiffen with the ladder).  GATED on the ENERGY CALIBRATION: when the reference grid already
    // sits at/above RelCutoffSafety() x kRungGateC x alpha_max the sharpest pair's tails are already inside
    // the calibrated budget and the rung would buy sub-mHa for 1.6-4x runtime (measured 2026-07-16 at the
    // Si anchors).  kRungGateC is a FIXED calibration constant, deliberately DECOUPLED from cutoffFactor:
    // the rung serves the pair->level rule (requirement max = RelCutoffSafety()*efine, no C in it) and the
    // absolute-kappa local-PP sweep -- when the (f2) raw-XC feed dropped the density-floor default C 8->2,
    // a gate borrowing C silently lost the rung at explicit Ecuts the kappa gate needs
    // (GPW.LocalPPKappaSelfConverged: analytic-vs-grid short 0.087 without the Ecut=20 rung at Ecut=10).
    static constexpr double kRungGateC=8.0;   // the C=8-era measured calibration, now standalone
    if (efine < itsLat->RelCutoffSafety()*kRungGateC*amax)
        levels.push_back(std::make_shared<const PW_Grid_Evaluator>(
                                grid->Recip(), rvec3_t(0,0,0), itsLat->RelCutoffSafety()*efine, itsRaster));
    for (const auto& g : levels)
    {
        levelN.push_back(g->FFTGrid());
        levelEcut.push_back(g->Ecut());
    }
}

// GRID DIAGNOSTIC (doc/GPWPlan §0e): one line per stored grid, printed at run start so GPW's grids can be
// lined up against CP2K's &MGRID log (NGRIDS/CUTOFF/REL_CUTOFF + per-level N).  |G|min is the smallest
// NONZERO |G| (the cell scale 2*pi/a); |G|max the ball edge (sqrt(2*Ecut)).
std::ostream& GPW_Evaluator::ReportGrid(std::ostream& os, const std::string& tag, const PW_Grid_Evaluator& g)
{
    double gmin=-1.0, gmax=0.0;
    for (const ivec3_t& m : g.Gs())
    {
        if (m.x==0 && m.y==0 && m.z==0) continue;
        const double gn=norm(g.GetGCartesian(m));
        if (gmin<0.0 || gn<gmin) gmin=gn;
        if (gn>gmax) gmax=gn;
    }
    const ivec3_t N=g.FFTGrid();
    return os<<"[GPW grid] "<<tag<<": N=("<<N.x<<","<<N.y<<","<<N.z<<") Ecut="<<g.Ecut()
             <<" Ha nG="<<g.size()<<" |G|min="<<gmin<<" |G|max="<<gmax<<"\n";
}

void GPW_Evaluator::ReportGrids(std::ostream& os) const
{
    os<<"[GPW basis] n="<<itsN<<" alpha_min="<<itsLat->MinExponent()<<" alpha_max="<<itsLat->MaxExponent()
      <<" cutoffFactor="<<itsCutoffFactor
      <<" raster="<<(itsRaster==RasterPolicy::BallOnly?"BallOnly":"AliasFree")
      <<" routing="<<(itsRelFieldSharp>=0.0?"HartreeOnly beta=":"HartreeXC beta=2/3*alpha_max=")
      <<(itsRelFieldSharp>=0.0?itsRelFieldSharp:(2.0/3.0)*itsLat->MaxExponent())<<" Ha-exp"
      <<" (auto density floor Ecut="<<itsCutoffFactor*itsLat->MaxExponent()<<" Ha)"<<std::endl;
    if (!itsFFT_R_G_Grids)
    {
        os<<"[GPW grid] DFT tier OFF (densityEcut=0): no grids"<<std::endl;
        return;
    }
    ReportGrid(os, "FFT rho<->G (density/collocation)", *itsFFT_R_G_Grids);
    EnsureLevels();
    // T1 {G}-STAR FOLD VISIBILITY (user 2026-08-15: "no hint anywhere of folding").  Under imposed ops
    // each level's ball partitions into stars; report the count + reduction beside nG.  HONESTY NOTE:
    // today the fold is CONSUMED only by the static local-PP sweeps (MakeLocalPP/MakeLocalPPLong T1
    // reduced evaluation); the per-iteration G-space consumers run unfolded -- doc/OpenWork.md item 5.
    std::vector<Symmetry::Lattice_3D::SymOp> sops;
    for (const auto& op : RecipSymOps()) sops.push_back({op.U, op.tau});
    for (size_t L=0;L<itsLevels.size();L++)
    {
        std::string tag="ladder L="+std::to_string(L);
        if (!sops.empty())
        {
            const size_t stars=Symmetry::Lattice_3D::FoldGVectors(itsLevels[L]->Gs(), sops).Reps();
            tag+=" stars="+std::to_string(stars);
        }
        if (L==0)                  tag+=" (== FFT grid; resolution reference)";
        else if (L>=itsNBaseLevels) tag+=" (top completion rung; local-PP uses base sub-ladder only)";
        ReportGrid(os, tag, *itsLevels[L]);
    }
    os<<"[GPW grid] local-PP integration: FULL ladder L=0.."<<itsLevels.size()-1
      <<" absolute REL_CUTOFF kappa="<<LocalPPRelCutoff()<<" Ha (e^{-kappa/2} pair tails)"<<std::endl;
    if (!sops.empty())
        os<<"[GPW grid] T1 {G}-star fold: "<<sops.size()<<" imposed ops; ACTIVE on the static local-PP"
            " sweeps (form factor at star reps); per-iteration G fields UNFOLDED (OpenWork item 5)"<<std::endl;
    os.flush();
}

void GPW_Evaluator::EmitGridsReport() const
{
    namespace rpt = qchem::report;
    rpt::json g;
    g["densityEcut"]  = itsCutoffFactor * itsLat->MaxExponent();   // the auto density-grid floor Ecut
    g["cutoffFactor"] = itsCutoffFactor;
    g["raster"]       = (itsRaster == RasterPolicy::BallOnly ? "BallOnly" : "AliasFree");
    g["routing"]      = (itsRelFieldSharp>=0.0 ? "HartreeOnly" : "HartreeXC");   // field-sharpness routing policy
    g["routingBeta"]  = (itsRelFieldSharp>=0.0 ? itsRelFieldSharp : (2.0/3.0)*itsLat->MaxExponent());
    if (itsFFT_R_G_Grids)                                          // null == DFT tier off (no grids)
    {
        EnsureLevels();
        // The {G}-star column (see ReportGrids): star counts per level under the imposed ops, so the
        // fold is VISIBLE in the grids table (nG vs Gstars = the T1 reduction).  Empty ops = no column.
        std::vector<Symmetry::Lattice_3D::SymOp> sops;
        for (const auto& op : RecipSymOps()) sops.push_back({op.U, op.tau});
        rpt::json ladder = rpt::json::array();
        for (size_t L = 0; L < itsLevels.size(); ++L)
        {
            const PW_Grid_Evaluator& lv = *itsLevels[L];
            const ivec3_t N = lv.FFTGrid();
            // role: L==0 is the density/collocation reference (== FFT grid); the top rung completes the
            // ladder above the local-PP sub-ladder; everything between is a coarser multigrid level.
            const char* role = (L == 0) ? "reference" : (L >= itsNBaseLevels ? "rung" : "coarse");
            rpt::json row;
            row["level"] = (long)L;
            row["N"]     = { N.x, N.y, N.z };
            row["ecut"]  = lv.Ecut();
            row["nG"]    = (long)lv.size();
            if (!sops.empty())
                row["Gstars"] = (long)Symmetry::Lattice_3D::FoldGVectors(lv.Gs(), sops).Reps();
            row["role"]  = role;
            ladder.push_back(row);
        }
        g["ladder"]  = ladder;
        if (!sops.empty())
            g["GstarFold"] = { { "ops", (long)sops.size() },
                               { "activeOn", "static local-PP sweeps (T1); per-iteration G fields unfolded" } };
        g["localPP"] = { { "kappa", LocalPPRelCutoff() } };
    }
    rpt::EmitSection("grids", g);
    // The lattice-sum economy readout (grids.latticeSums + the [lattice sums] console line): the numbers
    // that JUMP when diffuse functions are added -- the owner of the screens (the molecular evaluator)
    // reports them, so the eps values never leak through the face.
    itsLat->EmitLatticeSumReport(itsCell);
}

// Bloch sum of the Gaussian orbitals, chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R), over the COLLOCATION set (the
// orbital reach; == the overlap set unless collRcut decoupled it).  At Gamma every phase is 1, so the
// imaginary part is exactly zero (the sum reduces to the real molecular sum).
cvec_t GPW_Evaluator::Eval(const rvec3_t& r) const
{
    const double rr=itsCellRad+itsMaxReach;
    cvec_t v(itsN, dcmplx(0.0));
    for (size_t k=0;k<itsRc.size();k++)
    {
        const rvec3_t d=r-itsRc[k]-itsCellCtr;
        if (d.x*d.x+d.y*d.y+d.z*d.z > rr*rr) continue;   // image cannot reach r (every centre is in the cell)
        rvec_t chi=(*itsOrb)(r-itsRc[k]);
        for (size_t i=0;i<itsN;i++) v[i]+=itsPhaseC[k]*chi[i];
    }
    return v;
}

// Eval on a POINT SET.  The image sum is DELEGATED to the molecular seam (see the declaration): the
// enumeration is the seam's own (same eps-derived reach and radius rule this class used, so the
// untransformed table is bit-identical), and a transformed basis gets to transform once per point.
// Batched rather than per-point because the seam re-derives its offset list per call -- amortised over a
// block of points, paid per point if this were called pointwise.
mat_t<dcmplx> GPW_Evaluator::EvalMany(const rvec3vec_t& pts) const
{
    if (itsHomeOnly)     // molecule-in-a-box: image-free by definition, the raw home orbital widened
    {
        mat_t<dcmplx> P(pts.size(), itsN);
        for (size_t q=0;q<pts.size();q++)
        {
            const rvec_t chi=(*itsOrb)(pts[q]);
            for (size_t i=0;i<itsN;i++) P(q,i)=dcmplx(chi[i],0.0);
        }
        return P;
    }
    mat_t<dcmplx> Phi;
    itsLat->BlochPointValues(pts, CellPhase(), itsCell, Phi);
    return Phi;
}

cvec3vec_t GPW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    const double rr=itsCellRad+itsMaxReach;
    cvec3vec_t v(itsN, vec3_t<dcmplx>(dcmplx(0.0),dcmplx(0.0),dcmplx(0.0)));
    for (size_t k=0;k<itsRc.size();k++)
    {
        const rvec3_t d=r-itsRc[k]-itsCellCtr;
        if (d.x*d.x+d.y*d.y+d.z*d.z > rr*rr) continue;   // image cannot reach r
        rvec3vec_t g=itsOrb->Gradient(r-itsRc[k]);
        for (size_t i=0;i<itsN;i++)
            v[i]=v[i]+vec3_t<dcmplx>(itsPhaseC[k]*g[i].x, itsPhaseC[k]*g[i].y, itsPhaseC[k]*g[i].z);
    }
    return v;
}

// The periodic 1E matrices: analytic lattice sums, eps-CONVERGED inside the molecular seam (per-shell-pair
// internal enumeration -- the SAME screen the density collocation uses, so the two are ONE scheme by
// construction: the screened-complete Bloch operator).  S is PSD to eps << lambda_min (the screen drops
// only sub-eps terms of the exact PSD Gram).  KineticMatrix is <p^2> (no 1/2 -- the Hamiltonian applies
// it).  The home-only MODE returns the finite molecule's own (cached) matrices, widened to complex --
// the molecule-in-a-box configuration, by definition image-free.
namespace
{
template <class M> chmat_t Widen(const M& m)   // real symmetric -> complex Hermitian (zero imaginary part)
{
    const size_t n=m.rows();
    chmat_t c(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) c(i,j)=dcmplx(m(i,j),0.0);
    return c;
}
} //anon
chmat_t GPW_Evaluator::OverlapMatrix()                 const
{   qchem::report::Timed t("setup: analytic 1E lattice sums (S,T,V)");
    return itsHomeOnly ? Widen(itsOrb->Overlap())    : itsLat->MakeOverlap(CellPhase(),itsCell); }
chmat_t GPW_Evaluator::KineticMatrix()                 const
{   qchem::report::Timed t("setup: analytic 1E lattice sums (S,T,V)");
    return itsHomeOnly ? Widen(itsOrb->Kinetic())    : itsLat->MakeKinetic(CellPhase(),itsCell); }
chmat_t GPW_Evaluator::NuclearMatrix(const Structure* cl) const
{   qchem::report::Timed t("setup: analytic 1E lattice sums (S,T,V)");
    return itsHomeOnly ? Widen(itsOrb->Nuclear(cl))  : itsLat->MakeNuclear(CellPhase(),itsCell,cl); }

// The PP-quadrature integration mesh: a uniform lattice mesh whose Nyquist resolution follows the density
// cutoff (CreateIntegrationMesh's "GPW / Nyquist path").  The DFT tier (density grid) must be on: PP assembly
// only ever runs inside an SCF, which needs the density collocation grid anyway.
qcMesh::MeshParams GPW_Evaluator::PPMeshParams() const
{
    assert(itsFFT_R_G_Grids && "GPW_Evaluator: the external PP requires the DFT density grid (densityEcut!=0: <0 auto, >0 explicit)");
    qcMesh::MeshParams mp;
    mp.eCut=itsFFT_R_G_Grids->Ecut();
    return mp;
}

// The reciprocal {U|τ} face of the imposed ops (U = Wᵀ, the SymmetrizeGMap scatter convention) -- what the
// T1 {G}-star reduced structure-factor sweeps fold under.  Empty in, empty out (free run = plain sweep).
std::vector<Symmetry::Lattice_3D::ReciprocalOp> GPW_Evaluator::RecipSymOps() const
{
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> rops;
    rops.reserve(itsSymOps.size());
    for (const auto& op : itsSymOps) rops.push_back({Transpose(op.W), op.tau});
    return rops;
}

// Local PP: assembled in G-SPACE from the analytic form factor, IDENTICALLY to PW_Evaluator::LocalPotential-
// Matrix -- Vtilde(dG) = (1/Omega) Sum_a v_loc(Z_a,|dG|^2) e^{-i dG.tau_a}, dG=0 DROPPED, then <chi|Vtilde|chi>
// via OverlapMatrix (the collocation adjoint reconstructs V_loc(r) on the density grid and quadratures it).
// This inherits the PW G=0 / FormFactorG0-alignment convention EXACTLY, so the energy is box-independent (a
// real-space quadrature of the raw -Zion/r-tailed V_loc has a cell-size-dependent mean -> box drift).  The KB
// nonlocal (localized, no Coulomb tail, no G=0 issue) stays real-space (MakeSeparablePP).
chmat_t GPW_Evaluator::MakeLocalPP(const Structure* cl, const Pseudopotential::LocalPotential& loc, LocalPart part) const
{
    // In production this sweep runs for the LONG part only (PW_Hartree's V_long block; the short part is
    // analytic, MakeLocalPPShort below) -- the bucket name says so.  The full/short grid forms only run
    // for non-Gaussian PP models or the finite==lattice gates.
    qchem::report::Timed timed("setup: local-PP LONG integrate-back (grid sweep)");
    assert(itsFFT_R_G_Grids && "GPW_Evaluator: the local PP needs the density grid (densityEcut!=0: <0 auto, >0 explicit)");
    const UnitCell& B=itsFFT_R_G_Grids->Recip().GetCell();
    // ANALYTIC integrate-back of the G-space local-PP form factor, on the FULL ladder with the ABSOLUTE
    // pair->level rule (kappa, below).  V_loc is spectrally BROAD (unlike the smooth V_H/V_xc): at the
    // smooth-field calibration the mid pairs' e^{-2.5} tails against the deep well were a Ha-scale term
    // (Si SR Gamma -10.72 vs CP2K -7.115 with charge exact).  Fine-only (every pair on the fine grid) is
    // exact but UNTENABLE at a large fine grid (NaF's auto-160: an ultra-diffuse pair's box on
    // N=64 x ~180 screened offsets stalls the setup for hours) -- and unnecessary: an ultra-diffuse pair's OWN
    // spectrum kills the field tail, so it still belongs on a deep coarse level.  STATIC (framework-cached), so
    // this sweep is not a per-iteration cost; for Si SR (ladder {20,5}) scale 6 puts every pair on the fine
    // level == the validated fine-only numbers.
    auto Vt=[&](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0);        // drop dG=0 (alignment carries it)
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));
        double  g2=dG*dG;
        dcmplx  acc(0.0);                                            // form factor x structure factor
        for (Atom* a : *cl)                                          // full V_loc, or its long / short piece
        {
            double f = part==LocalPart::Long  ? loc.FormFactorLong (a->itsZ,g2)
                     : part==LocalPart::Short ? loc.FormFactorShort(a->itsZ,g2)
                     :                          loc.FormFactor     (a->itsZ,g2);
            acc += f*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        }
        return acc/itsFFT_R_G_Grids->Volume();
    };
    EnsureLevels();
    // FULL ladder + the ABSOLUTE assignment rule (doc/GPWPlan.md 0e-PP step (a), 2026-07-22).  The old
    // relCutoffScale (a multiplier on the RELATIVE smooth-field rule) could not make this sweep
    // standalone-exact: the relative rule's tail bound scales with ecut_ref/alpha_max, so the split long /
    // short pieces each carried ~0.5 Ha band-limit errors that only CANCELLED in the smooth sum -- fatal
    // once the short goes analytic (step (b)).  The absolute rule req = kappa*(alpha_i+alpha_j) bounds
    // EVERY pair's spectral tail by e^{-kappa/2} uniformly (kappa=30 Ha -> e^{-15}), independent of the
    // field's sharpness -- CP2K's gaussian_gridlevel REL_CUTOFF (60 Ry), and the reason its numeric-but-
    // smooth V_long integration is sub-mHa.  The TOP COMPLETION RUNG is included (full ladder): the
    // sharpest pairs' requirement kappa*2*alpha_max wants it; pairs beyond the finest level fall back to
    // the finest present (CP2K's gridlevel=1 default).  GPW_LOCALPP_RELCUTOFF overrides kappa (the
    // self-convergence verification: kappa=60 -> e^{-30} must match kappa=30 to tolerance).
    const double kappa=LocalPPRelCutoff();
    const size_t K = itsLevels.size();
    // T1 {G}-star fold (doc/SymmetryUpgradePlan.md §6): under the §3-imposed ops the static V_loc field is
    // exactly symmetric (the ops were detected from these very atoms), so the form-factor sum runs at star
    // representatives only, members filled by the glide identity.  Empty ops (free run) = the plain sweep.
    const auto rops=RecipSymOps();
    std::vector<rvec_t> V_L(K);
    size_t gRaw=0, gReps=0;                        // Step 0b: the achieved {G}-star reduction, summed over levels
    for (size_t L=0;L<K;L++)
    {
        size_t nr=0;
        ΔG_Map vmapL = EvaluateSymmetricGMap(itsLevels[L]->Gs(), rops, Vt, &nr);   // restrict to level L's {G}
        gRaw+=itsLevels[L]->Gs().size(); gReps+=nr;
        V_L[L]=itsLevels[L]->RhoOnGrid(vmapL);
    }
    qchem::report::EmitFold("V_loc {G}-star", rops.size(), gRaw, gReps);
    // THE COLLOCATED FIELD's OWN SHARPNESS (doc/GPWPlan1.md 4b root-cause fix).  V_loc (long OR short) decays in
    // G as e^{-G^2 rloc^2/2}, so its effective Gaussian exponent is beta = 1/(2 rloc^2) -- exactly the short-
    // range Gaussian's alpha.  The integrand chi_i*V*chi_j is a product of Gaussians (exponents ADD), so the
    // pair->level rule must resolve alpha_i+alpha_j+beta, NOT the pair alone: without beta a DIFFUSE pair against
    // a SHARP PP well (small rloc, e.g. F) lands on a coarse level that cannot resolve the well -> a spurious
    // low diffuse eigenvalue (the NaF diving-ghost).  beta = max over species = the SHARPEST (smallest rloc) PP,
    // since V_loc sums over all atoms.  A soft PP (large rloc) or a non-Gaussian model -> beta=0 = old behaviour.
    double beta=0.0;
    if (const auto* gauss=dynamic_cast<const Pseudopotential::LocalPotential_Gaussian*>(&loc))
        for (Atom* a : *cl)
        {
            const auto terms=gauss->ShortRangeGaussian(a->itsZ);
            if (!terms.empty()) beta=std::max(beta, terms[0].alpha);   // alpha = 1/(2 rloc^2)
        }
    const bool timeIt=(std::getenv("GPW_LOCALPP_RELCUTOFF")!=nullptr);
    auto t0=std::chrono::steady_clock::now();
    chmat_t h=itsLat->IntegratePotential(V_L, CellPhase(), itsCell, itsLevelN, itsLevelEcut, kappa, nullptr, beta);
    if (timeIt)
    {
        double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t0).count();
        std::cerr<<"[MakeLocalPP part="<<int(part)<<" kappa="<<kappa<<" K="<<K
                 <<"] IntegratePotential "<<ms<<" ms"<<std::endl;
    }
    return h;
}

// The LONG-range (softened-Coulomb) local-PP matrix, the CP2K split's Hartree-folded piece (doc/GPWPlan.md
// 0e-PP / 0i increment 3 -- the CUSTOM V_loc G-BALL, 2026-08-01).  V_long is a STATIC EXTERNAL field of
// effective Gaussian exponent beta = 1/(2 rloc^2): unlike the per-iteration KS fields -- whose per-level
// band-limit is EXACT by adjoint-consistency with the density collocation -- its level restriction is a
// plain truncation error, so it is ENTITLED TO ITS OWN G-BALL sized by the integrand's spectra (the user's
// additive principle vlocEcut ~ 2 alpha_max + beta, sharpened to the exact HARMONIC truncation bound
// 2 ln(1/eps) p beta/(p+beta) per pair -- StaticFieldPairLevels, doc there).  The assembly: the block's own
// ladder plus ONE custom top level at the sharpest pair's requirement, each pair routed by the harmonic
// rule.  This is BOTH the accuracy fix and the speed fix over the retired kappa sweep:
//  - the smooth-fold experiment (0i increment 2, 2026-07-31) put mid pairs on the Ecut-80 ball where the
//    field tail e^{-g2 rloc^2/2} ~ e^{-3.8} cost 4.6 mHa on the diffuse NaF; the harmonic rule routes them
//    to the custom top ball (e^{-2 ln(1/eps)} class);
//  - the kappa sweep's max(p,beta) requirement mis-routed the DIFFUSE pairs onto the completion rung
//    (req = kappa*beta ~ 315 Ha even for p~0.2) -- huge boxes x ~791 offsets = the 180 s dominant
//    diffuse-NaF setup bucket; under the harmonic rule their own spectrum kills the field tail and they
//    fall to deep coarse levels (tiny work), so the on-the-fly sweep costs seconds.
// TWO invariants: screenD stays NULL -- a STATIC (density-independent) block must not freeze one
// iteration's D-aware active set -- and screenD==null makes the phase-independent B_ij memo apply, so
// multi-k pays the assembly once.  eps defaults 1e-5 (the kappa-sweep-parity class: its rung-160 tails
// were e^{-8.5}, measured sub-mHa vs CP2K); GPW_VLOC_EPS overrides for the self-convergence check.
// GPW_LONG_SWEEP=1 = the kappa-sweep path (A/B verification instrument); non-Gaussian local models (no
// closed beta) also fall back to it.
chmat_t GPW_Evaluator::MakeLocalPPLong(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    // beta = the long part's effective exponent = the SHARPEST species' 1/(2 rloc^2) (V_long sums over atoms)
    double beta=0.0;
    if (const auto* gauss=dynamic_cast<const Pseudopotential::LocalPotential_Gaussian*>(&loc))
        for (Atom* a : *cl)
        {
            const auto terms=gauss->ShortRangeGaussian(a->itsZ);
            if (!terms.empty()) beta=std::max(beta, terms[0].alpha);   // alpha = 1/(2 rloc^2)
        }
    static const bool oldSweep=(std::getenv("GPW_LONG_SWEEP")!=nullptr);
    if (oldSweep || beta<=0.0) return MakeLocalPP(cl, loc, LocalPart::Long);
    qchem::report::Timed timed("setup: local-PP LONG (custom G-ball)");
    assert(itsFFT_R_G_Grids && "GPW_Evaluator: the local PP needs the density grid (densityEcut!=0)");
    EnsureLevels();
    static const double eps=[]{const char* s=std::getenv("GPW_VLOC_EPS"); return s?std::atof(s):1e-5;}();
    const double lnE=-std::log(eps);
    // The custom ladder: the block's levels + ONE top level at the sharpest pair's harmonic requirement
    // (2 alpha_max against beta) -- appended LAST so ecut_L[0] stays the resolution reference.  Skipped when
    // an existing level (the completion rung) already satisfies it.  The requirement SATURATES at
    // 2 lnE beta for a hard basis (the field's own bandwidth is the ceiling), so the ball never runs away.
    const double p2=2.0*itsLat->MaxExponent();
    const double ecutTop=2.0*lnE*p2*beta/(p2+beta);
    std::vector<ivec3_t> N_L=itsLevelN;
    std::vector<double>  ecut_L=itsLevelEcut;
    std::vector<std::shared_ptr<const PW_Grid_Evaluator>> levels=itsLevels;
    double emax=0.0;
    for (double e : ecut_L) emax=std::max(emax,e);
    if (ecutTop>emax)
    {
        auto g=std::make_shared<const PW_Grid_Evaluator>(itsFFT_R_G_Grids->Recip(), rvec3_t(0,0,0), ecutTop, itsRaster);
        levels.push_back(g); N_L.push_back(g->FFTGrid()); ecut_L.push_back(g->Ecut());
    }
    const std::vector<size_t> lv=itsLat->StaticFieldPairLevels(ecut_L, beta, lnE);
    // Vtilde(dm) = (1/Omega) Sum_a f_long(Z_a,|dG|^2) e^{-i dG.tau_a}, dG=0 DROPPED (the E_alphaZ
    // alignment carries the mean -- the same convention as MakeLocalPP's form-factor closure).
    const UnitCell& B=itsFFT_R_G_Grids->Recip().GetCell();
    auto Vt=[&](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0);
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));
        double  g2=dG*dG;
        dcmplx  acc(0.0);                                            // form factor x structure factor
        for (Atom* a : *cl) acc += loc.FormFactorLong(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        return acc/itsFFT_R_G_Grids->Volume();
    };
    const size_t K=levels.size();
    const auto rops=RecipSymOps();                                    // T1 {G}-star fold, as in MakeLocalPP
    std::vector<rvec_t> V_L(K);
    size_t gRaw=0, gReps=0;
    for (size_t L=0;L<K;L++)
    {
        size_t nr=0;
        ΔG_Map vmapL = EvaluateSymmetricGMap(levels[L]->Gs(), rops, Vt, &nr);   // level L's {G} (low-pass)
        gRaw+=levels[L]->Gs().size(); gReps+=nr;
        V_L[L]=levels[L]->RhoOnGrid(vmapL);
    }
    qchem::report::EmitFold("V_loc-long {G}-star", rops.size(), gRaw, gReps);
    return itsLat->IntegratePotential(V_L, CellPhase(), itsCell, N_L, ecut_L,
                                      0.0, nullptr, 0.0, itsRelFieldSharp, &lv);   // explicit harmonic routing; NO D screen
}

// The SHORT-range local-PP matrix, ANALYTIC via the closed Gaussian form (increment 2, doc/GPWPlan.md 0e-PP):
// V_short(r) = Sum_t c_t r^{2n_t} e^{-alpha r^2} at each atom, so <chi_i|V_short|chi_j> is a lattice-summed
// 3-centre Overlap3C (LatticeSum1E::MakeLocalGaussian) -- NO grid sweep (the compact core needs no
// relCutoffScale stiffening / giant boxes).  Falls back to the grid sweep for a model whose short part is not
// exposed in closed Gaussian form.  The per-species operator g_Z = Sum_t c_t r^{2n_t} e^{-alpha r^2} is a
// Cartesian-Gaussian function: r^{2n} = MultiplyR2(1,n) (the l=0 s-harmonic is the constant), all sharing the
// one exponent alpha = 1/(2 rloc^2); each atom places g_Z at its centre (the OPERATOR slot of Overlap3C).
chmat_t GPW_Evaluator::MakeLocalPPShort(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    const auto* gauss=dynamic_cast<const Pseudopotential::LocalPotential_Gaussian*>(&loc);
    if (!gauss) return MakeLocalPP(cl, loc, LocalPart::Short);   // grid fallback (non-Gaussian short part)
    qchem::report::Timed timed("setup: local-PP SHORT (analytic lattice sum)");
    auto opForZ=[gauss](int Z)->Molecule::LatticeSum1E::GaussianFunction
    {
        Molecule::LatticeSum1E::GaussianFunction g;
        const auto terms=gauss->ShortRangeGaussian(Z);
        g.alpha = terms.empty() ? 1.0 : terms[0].alpha;
        for (const auto& rt : terms)
        {
            assert(rt.alpha==g.alpha && "GPW short PP: a species' Gaussian terms must share one exponent (1/2 rloc^2)");
            auto mono=MultiplyR2({ qchem::Math::CartTerm{ qchem::Math::Monomial{0,0,0}, 1.0 } }, rt.n);
            for (auto& t : mono) { t.c*=rt.c; g.terms.push_back(t); }
        }
        return g;
    };
    chmat_t V = itsHomeOnly ? itsLat->MakeLocalGaussian(cl, opForZ)
                            : itsLat->MakeLocalGaussian(CellPhase(), itsCell, cl, opForZ);
    // PERIODIC G=0 CONVENTION (the 5.7% wiring bug, 2026-07-22): the grid sweep's form-factor closure DROPS
    // dG=0 -- the cell mean of V_short is carried by the E_alphaZ term (PW_Pseudo, FormFactorG0Short), not the
    // matrix.  The analytic lattice sum computes the FULL integral including that mean, so subtract it
    // (V-bar * S, with S this block's own Bloch overlap) or the mean would be DOUBLE-counted against E_alphaZ.
    // A genuinely finite Structure keeps the full integral (no neutralising background, no E_alphaZ) -- the
    // same isFinite() physics decision PW_Pseudo makes.
    if (!cl->isFinite())
    {
        double vbar=0.0;
        for (Atom* a : *cl) vbar += loc.FormFactorG0Short(a->itsZ);
        vbar /= itsCell.GetCellVolume();
        V -= vbar*OverlapMatrix();
    }
    return V;
}

// KB separable nonlocal PP: accumulate the rank-1 Hermitian D|b><b| over atoms/projectors/m, with the BLOCH
// projection vector b_i^k = <chi_i^k | beta_p Y_lm> = Sum_R e^{ik.R} integral chi_i^k(r)* beta(r-(R_a-R)) d3r,
// mesh-quadratured over the home cell.
//
// THE BRA IS THE PERIODIC (BLOCH-SUMMED) ORBITAL chi_i^k = Eval, NOT the raw home orbital *itsOrb (fixed
// 2026-07-09, doc/GPWPlan.md TODO 1a).  A lattice-site atom (the FCC corner atom at 0) has its Gaussian split
// by the cell boundary; the RAW orbital on the single-cell integration mesh keeps only the in-cell fraction,
// so <chi_i|beta> lost the wrapped tail and the KB energy was translation-variant by ~16 Ha (images did NOT
// cure it -- summing the PROJECTOR images cannot restore the ORBITAL's missing tail).  The Bloch orbital
// brings the wrapped tail back onto the cell mesh via chi_i's neighbouring-cell image, so integral_cell
// chi_i^k* beta_a^k = the correct all-space projection (consistent with the analytic all-space overlap S).
// Needs the orbital reach (Rcut/collRcut>0) large enough to wrap; at Rcut=0 Eval==raw orbital (no wrap) so the
// committed Gamma/atom-in-box anchors, all Rcut=0, are UNCHANGED.  Per-term gate: CP2K Nonlocal PP = +0.9406 Ha.
//
// chi_i^k is precomputed on the mesh ONCE (Eval re-sums the image set per point, so evaluating it inside the
// per-projector-image loop would be O(images^2)); the projection then reuses it.  V_ij = Sum D b_i conj(b_j)
// is Hermitian by construction.
chmat_t GPW_Evaluator::MakeSeparablePP(const Structure* cl, const Pseudopotential::SeparablePotential_R& sep) const
{
    // The lumped matrix == the exact sum of the per-channel decomposition (ONE assembly, two faces).
    chmat_t H=blazem::zeroH<dcmplx>(itsN);
    for (const auto& lH : MakeSeparablePPByL(cl,sep)) H+=lH.second;
    return H;
}

std::map<int,chmat_t> GPW_Evaluator::MakeSeparablePPByL(const Structure* cl, const Pseudopotential::SeparablePotential_R& sep) const
{
    qchem::report::Timed timed("setup: separable PP (KB)");
    const size_t n=itsN;
    std::map<int,mat_t<dcmplx>> Vl;                  // one raw accumulator per projector channel l
    auto Vof=[&Vl,n](int l)->mat_t<dcmplx>&
    {
        auto it=Vl.find(l);
        if (it==Vl.end()) it=Vl.emplace(l,mat_t<dcmplx>(n,n,dcmplx(0.0))).first;
        return it->second;
    };
    auto pack=[n](std::map<int,mat_t<dcmplx>>& vl)   // project each channel to Hermitian (upper; real diag)
    {
        std::map<int,chmat_t> out;
        for (auto& lV : vl)
        {
            const mat_t<dcmplx>& V=lV.second;
            chmat_t H=blazem::zeroH<dcmplx>(n);
            for (size_t i=0;i<n;i++)
            {
                H(i,i)=dcmplx(std::real(V(i,i)),0.0);
                for (size_t j=i+1;j<n;j++) H(i,j)=V(i,j);
            }
            out.emplace(lV.first,std::move(H));
        }
        return out;
    };
    // ANALYTIC path (2026-07-15): a GTH/HGH projector is polynomial x Gaussian, so when the model exposes its
    // closed Gaussian form (SeparablePotential_Gaussian) the Bloch projection is an ANALYTIC lattice-summed
    // overlap -- b_i = Sum_R e^{-ik.R} <chi_i | beta Y_lm at tau_a - R> maps onto the molecular seam's
    // Sum_r phases[r] <chi_i | g(.-Rs[r])> with Rs = -itsRc, phases = conj(itsPhaseC) (the SAME conjugated
    // image phase as the mesh path; the all-space tiling derivation in the block comment below).  beta Y_lm =
    // Sum_t c_t r^{2n_t} e^{-alpha_t r^2} x [r^l Y_lm] with r^l Y_lm a degree-l Cartesian polynomial
    // (YlmCartesian, pinned to this file's RealYlm convention) -- one GaussianFunction per radial term.
    // EXACT (no mesh, no quadrature error) and ~O(n x images) 2-centre integrals instead of the mesh sweep
    // (NaF: the 358k-point Eval quadrature was >20 min of setup; this is milliseconds).  Models without the
    // closed-Gaussian face keep the mesh quadrature below.
    if (const auto* gsep=dynamic_cast<const Pseudopotential::SeparablePotential_Gaussian*>(&sep))
    {
        // b_i = <chi_i^k | beta at tau> = Sum_n phase(n) <chi_i | g(.-tau-R_n)> -- the seam enumerates the
        // series internally per (chi_i, g).  (The historical (-Rs, conj-phase) form was an artifact of the
        // pre-built inversion-symmetric list: substituting m=-n gives conj(phase(-m))=phase(m), i.e. the
        // PLAIN phase oracle with g imaged at +R_m.  Gate: AnalyticSeparablePPMatchesMesh == mesh.)
        // Home-only mode: the single finite term <chi_i|g> (the raw home orbital -- the box configuration).
        auto phase=CellPhase();
        for (size_t a=0; a<cl->GetNumAtoms(); a++)
        {
            const Atom* at=(*cl)[a];
            int Z=at->itsZ;
            for (size_t p=0; p<sep.NumProjectors(Z); p++)
            {
                int    l=sep.AngularMomentum(Z,p);
                double D=sep.Coefficient    (Z,p);
                auto   radial=gsep->BetaGaussian(Z,p);
                for (int m=-l; m<=l; m++)
                {
                    auto ylm=YlmCartesian(l,m);
                    cvec_t b(n, dcmplx(0.0));
                    for (const auto& rt : radial)               // one seam call per radial term c r^{2n} e^{-ar^2}
                    {
                        Molecule::LatticeSum1E::GaussianFunction g;
                        g.center=at->itsR;
                        g.alpha =rt.alpha;
                        g.terms =MultiplyR2(ylm, rt.n);
                        for (auto& t : g.terms) t.c*=rt.c;
                        cvec_t bt = itsHomeOnly ? itsLat->MakeOverlap(g) : itsLat->MakeOverlap(phase, itsCell, g);
                        for (size_t i=0;i<n;i++) b[i]+=bt[i];
                    }
                    mat_t<dcmplx>& V=Vof(l);
                    for (size_t i=0;i<n;i++)
                        for (size_t j=0;j<n;j++) V(i,j)+=D*b[i]*std::conj(b[j]);
                }
            }
        }
        return pack(Vl);
    }

    qcMesh::Mesh mesh=cl->CreateIntegrationMesh(PPMeshParams());
    const rvec3vec_t& R=mesh.Points();
    const rvec_t&     W=mesh.Weights();
    size_t npts=mesh.size();

    // The KB projector beta_p is COMPACTLY supported (localized at each atom), so at a large image reach most
    // images place the projector centre far outside the cell (contributing exp(-large)~0), and even a near
    // image touches the cell only within its support radius r_beta.  SCREEN both: (a) skip an image whose
    // projector centre cannot reach the mesh (bounding sphere), (b) skip a mesh point beyond r_beta.  The Bloch
    // orbital chi_i^k on the mesh (Eval, the dominant cost) is then needed ONLY at the surviving points, so it
    // is computed LAZILY.  Screening is at 1e-10 of the projector peak -> numerically exact (anchors unmoved),
    // and cuts this build ~5x (the 133-image Rcut=2a sum collapses to the ~dozen images that touch the cell).
    rvec3_t ctr(0,0,0); for (size_t k=0;k<npts;k++) ctr=ctr+R[k]; ctr=ctr/double(npts);
    double  rad=0.0;    for (size_t k=0;k<npts;k++) rad=std::max(rad, norm(R[k]-ctr));

    mat_t<dcmplx> Phi(npts,n,dcmplx(0.0));   // chi_i^k on the mesh; filled lazily at screened-in points only
    std::vector<char> havePhi(npts,0);
    auto PhiAt=[&](size_t k){ if (!havePhi[k]) { cvec_t v=Eval(R[k]); for (size_t i=0;i<n;i++) Phi(k,i)=v[i]; havePhi[k]=1; } };

    for (size_t a=0; a<cl->GetNumAtoms(); a++)
    {
        const Atom* at=(*cl)[a];
        int Z=at->itsZ;
        for (size_t p=0; p<sep.NumProjectors(Z); p++)
        {
            int    l=sep.AngularMomentum(Z,p);
            double D=sep.Coefficient    (Z,p);
            const double rBeta=BetaSupportRadius(sep,Z,p);   // projector support radius (screening cutoff)
            for (int m=-l; m<=l; m++)
            {
                cvec_t b(n, dcmplx(0.0));            // b_i = Sum_R e^{-ik.R} integral_cell chi_i^k* beta(.-(R_a-R))
                for (size_t r=0; r<itsRc.size(); r++)   // orbital reach: the collocation image set
                {
                    rvec3_t c=at->itsR-itsRc[r];                 // this image's projector centre
                    if (norm(c-ctr) > rad+rBeta) continue;       // (a) image too far to touch the mesh
                    BetaYlmField beta(c,sep,Z,p,l,m);
                    // The projector-image phase is CONJUGATED (e^{-ik.R}), NOT e^{+ik.R}.  b_i = <chi_i^k|beta_home>
                    // = integral_allspace chi_i^k* beta(.-tau_a); tiling all-space into home cells (int_all f =
                    // Sum_R int_home f(.+R)) and applying the Bloch law chi_i^k(r+R)=e^{ik.R}chi_i^k(r) puts
                    // e^{-ik.R} on conj(chi_i^k) against the R-shifted projector beta(.-(tau_a-R)).  At Gamma /
                    // half-integer k every phase is +-1 (self-conjugate) so this was inert -- the FIRST genuinely
                    // complex k (e.g. k=1/4) is where +ik.R halved the KB trace (Vnl 42->22) and over-bound.
                    dcmplx ph=std::conj(itsPhaseC[r]);
                    for (size_t k=0;k<npts;k++)
                    {
                        if (norm(R[k]-c) > rBeta) continue;      // (b) mesh point beyond the projector support
                        double bw=beta(R[k])*W[k];
                        PhiAt(k);
                        for (size_t i=0;i<n;i++) b[i]+=ph*std::conj(Phi(k,i))*bw;
                    }
                }
                mat_t<dcmplx>& V=Vof(l);
                for (size_t i=0;i<n;i++)
                    for (size_t j=0;j<n;j++) V(i,j)+=D*b[i]*std::conj(b[j]);
            }
        }
    }
    return pack(Vl);
}

// Cache key: the molecular basis's geometry-aware ID pins the radials + centres (atom positions); the CELL
// VECTORS enter via volume + max edge (the atom positions alone cannot distinguish two cells); k + the
// home-only mode + the density-grid cutoff distinguish the periodic block.
std::string GPW_Evaluator::IDFragment() const
{
    // Include the density-grid cutoff: the collocation tensor (Repulsion3C/Overlap3C) is built on that grid, so
    // the framework cache (keyed by BasisSetID) must distinguish GPW bases that differ only in densityEcut.
    // The T3 stream-fold state (§6b/T3.2): a folded evaluator collocates REDUCED tensors (completed by the
    // caller-side group-average), so folded and unfolded runs over the same geometry must never share
    // framework-cached tensor closures (measured: the A/B gate's second run silently replayed the first
    // run's unfolded closures until this key landed).
    const size_t sfold = itsLat ? itsLat->StreamFoldOrder() : 0;
    return "|mol="+itsOrb->BasisSetID()
         +"|k="+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
         +"|cell="+std::to_string(itsCell.GetCellVolume())+","+std::to_string(itsCell.GetMaximumCellEdge())
         +(itsHomeOnly?"|home":"")
         +(sfold?"|sfold="+std::to_string(sfold):"")
         +"|dEcut="+(itsFFT_R_G_Grids?std::to_string(itsFFT_R_G_Grids->Ecut()):std::string("0"));
}

} //namespace
