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
module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.Blaze;       // rvec_t, rmat_t, rsmat_t, blazem::zeroH<dcmplx>
import qchem.Vector3D;    // vec3_t + rvec3_t / rvec3vec_t arithmetic (r - R, componentwise add)
import qchem.Mesh.Quadrature; // qcMesh::WeightedOverlap / Overlap (the real-space PP quadrature primitives)
import qchem.ScalarFunction;  // ScalarFunction<double> (the V_loc / beta*Ylm fields handed to the quadrature)
import qchem.Math;            // norm, Pi, sqrt (the real spherical harmonics for the KB projectors)
import qchem.Math.Angular;    // Monomial/CartTerm/SphericalShell (the analytic KB Cartesian expansion)
import qchem.Structure;       // Structure / Atom (the PP centres + Z, and CreateIntegrationMesh)
import qchem.UnitCell;        // UnitCell (the direct cell for CollocateDensity / IntegratePotential grid<->cell)

namespace qchem::BasisSet::Lattice_3D
{

namespace
{
// Build a lattice-translation set {R} (Cartesian, origin first) and its matching Bloch phases {e^{ik.R}} for
// a cutoff radius -- the {R}+{phase} weighted point set (future: one qcMesh cMesh).  phase = exp(2 pi i k.n)
// with n the integer cell index (convention-safe).  Rcut<=0 -> the home cell only (origin, phase 1).
void BuildImages(const UnitCell& cell, double Rcut, const rvec3_t& kFrac,
                 std::vector<rvec3_t>& R, cvec_t& phase)
{
    blazem::VecBuilder<dcmplx> ph;
    if (Rcut>0.0)
        for (const auto& n : cell.CellsInSphere(Rcut))
        {
            R.push_back(cell.ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z))));
            double kn=kFrac.x*n.x + kFrac.y*n.y + kFrac.z*n.z;         // k_frac . n
            ph.Append(std::exp(dcmplx(0.0, 2.0*Pi*kn)));              // e^{ik.R} = e^{2 pi i k_frac . n}
        }
    else
    {
        R.push_back(rvec3_t(0,0,0));
        ph.Append(dcmplx(1.0,0.0));
    }
    phase=ph.take();
    assert(!R.empty() && R.size()==phase.size());
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
                             double densityEcut, const rvec3_t& kFrac, double Rcut, double collRcut,
                             double cutoffFactor)
    : itsMol(std::move(mol))
    , itsk(kFrac)
    , itsCell(cell)
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

    // The 1E/orbital image set {R}+{e^{ik.R}}.  THREE modes (mirroring densityEcut):
    //   Rcut < 0  AUTO (recommended): the enumeration radius is DERIVED from the basis -- the magnitude screen
    //             keeps a pair (i, j@R) only within reach_i+reach_j <= 2 sqrt(-ln eps/alpha_min), and the two
    //             centres can sit up to a cell span apart, so 2*reach + 2*maxEdge enumerates EVERYTHING the
    //             screen can keep.  Screening then prunes it sparse (cost ~ the surviving terms, not Rcut^3).
    //             This removes the user-facing Rcut knob: the analytic collocation already enumerates its own
    //             offsets (ForImageOffsets), and with AUTO the 1E sums / Eval / KB projector do too.
    //   Rcut = 0  home cell only (the FINITE-molecule mode: reproduces the finite matrices exactly).
    //   Rcut > 0  explicit enumeration radius (the legacy knob; screening still prunes inside it).
    // CellsInSphere is inversion-symmetric, so the lattice-sum matrices come out Hermitian (real at Gamma).
    // collRcut<=0 reuses the overlap set (backward-compatible: Gamma/finite runs collocate on the same set).
    double rcutEff=Rcut;
    if (Rcut<0.0)
    {
        const double reach=std::sqrt(-std::log(1e-10)/itsLat->MinExponent());  // == the kernels' kScreenEps reach
        rcutEff=2.0*reach + 2.0*cell.GetMaximumCellEdge();
    }
    BuildImages(cell, rcutEff, itsk, itsR, itsPhase);
    // Per-image screening data for Eval/EvalGradient: every orbital centre sits inside the cell, so the image
    // R contributes at a point r only if |(r-R) - cellCentre| <= cellRadius + maxReach.  One norm per image
    // (vs n orbital evaluations) -- without it the AUTO enumeration radius makes the Bloch sum O(|images|)
    // dense at every KB mesh point (the 1E pair radius is ~2x what the single-orbital sum needs).
    itsMaxReach=std::sqrt(-std::log(1e-10)/itsLat->MinExponent());
    itsCellCtr =cell.ToCartesian(rvec3_t(0.5,0.5,0.5));
    itsCellRad =0.0;
    for (int cx=0;cx<=1;cx++) for (int cy=0;cy<=1;cy++) for (int cz=0;cz<=1;cz++)
    {
        rvec3_t corner=cell.ToCartesian(rvec3_t(double(cx),double(cy),double(cz)));
        itsCellRad=std::max(itsCellRad, norm(corner-itsCellCtr));
    }
    if (collRcut>0.0) BuildImages(cell, collRcut, itsk, itsRc, itsPhaseC);
    else { itsRc=itsR; itsPhaseC=itsPhase; }

    // The DFT tier's density/collocation grid: GPW's ONLY grid cutoff (no orbital/wavefunction cutoff -- the
    // Gaussians are analytic).  IMPORTANT: densityEcut is a DENSITY-scale quantity; the sharpest feature is the
    // PRODUCT of the two tightest primitives (a Gaussian of exponent 2*alpha_max), so at a fixed tolerance it
    // needs ~2x the cutoff that resolves a single orbital exp(-alpha_max r^2).  cutoffFactor (default 4) FOLDS
    // that x2 into a tolerance calibrated on the DENSITY charge loss, not an orbital tail (F alpha_max=40:
    // Ecut=40 loses >5 e-, =120 holds 7.997, CP2K converges ~200), so 4*alpha_max lands in the good regime.
    //   densityEcut < 0 : AUTOMATIC (recommended) -- floor = cutoffFactor*alpha_max, from the basis (no user Ha).
    //   densityEcut = 0 : DFT tier OFF (1E-only; no grid).
    //   densityEcut > 0 : EXPLICIT -- honoured as given, but WARN on cerr if below the floor (the caller insisted
    //                     on an under-resolved grid: charge leaks off-grid -- we don't hide it, but we don't
    //                     silently override the explicit choice either).  For SIPP alpha_max=2 the floor is 8,
    //                     below every committed Si anchor's densityEcut (>=10), so those are unchanged.
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
        itsGrid=std::make_shared<const PW_Grid_Evaluator>(
                    ReciprocalLattice(cell.MakeReciprocalCell()), rvec3_t(0,0,0), ecut);
    }
}

// The Bloch phase closure the analytic kernels call back per screened cross-cell offset: e^{2 pi i k_frac.n}
// (n the INTEGER cell index -- convention-safe, same convention as BuildImages).  Gamma: the constant 1.
Molecule::LatticeSum1E::cellphase_t GPW_Evaluator::CellPhase() const
{
    const rvec3_t k=itsk;
    return [k](const ivec3_t& n)->dcmplx
    {
        double kn=k.x*n.x + k.y*n.y + k.z*n.z;
        return std::exp(dcmplx(0.0, 2.0*Pi*kn));
    };
}

// The matrix-free density->rho-tilde / ->V_H map (the G_ERI3::apply realization -- the CP2K analytic
// collocation, doc/GPWPlan.md S0).  The density is collocated ANALYTICALLY per (pair, screened cross-cell
// offset) on compact exp-tail boxes, modulo-wrapped, each pair on the coarsest REL_CUTOFF level resolving its
// product exponent (LatticeSum1E::CollocateDensity); each level is FFT'd and rho-tilde is combined NESTED in
// G-space (a coarse level contributes over its own {G} -- beyond it the pair's spectrum is below tolerance by
// construction).  rho-tilde(G)=FFT[rho]/Npts (the grid quadrature's 1/Omega and Omega/Npts cancel).  For
// Coulomb the diagonal Poisson kernel 4pi/|G|^2 is folded in (V_H = 4pi rho-tilde/G^2; G=0 -> 0).  The closure
// keeps the level grids + molecular basis alive (captured shared_ptrs) since it lives in the framework-cached
// G_ERI3.
std::function<ΔG_Map(const chmat_t&)> GPW_Evaluator::MakeCollocator(bool coulomb) const
{
    EnsureLevels();
    const UnitCell A = itsCell;                             // the direct cell (stored; see the member note)
    auto levels = itsLevels;                                // shared_ptrs -> the closure keeps the grids alive
    auto N_L    = itsLevelN;
    auto ecut_L = itsLevelEcut;
    auto mol    = itsMol;                                   // shared_ptr -> keep the molecular basis (lat) alive
    const Molecule::LatticeSum1E* lat = itsLat;
    auto phase  = CellPhase();
    ReciprocalLattice recip = itsGrid->Recip();
    if (!itsCollocMemo) itsCollocMemo=std::make_shared<CollocMemo>();
    auto memo   = itsCollocMemo;                            // ONE memo across the Coulomb + overlap closures
    return [A,levels,N_L,ecut_L,mol,lat,phase,recip,coulomb,memo](const chmat_t& D) -> ΔG_Map
    {
        // Same D as the last collocation (the sibling tensor's call this iteration): replay the level
        // densities.  EXACT equality (bit-identical or recompute); phase/cell/ladder are fixed per evaluator.
        auto sameD=[&]() -> bool
        {
            if (!memo->valid || memo->D.rows()!=D.rows()) return false;
            for (size_t i=0;i<D.rows();i++)
                for (size_t j=i;j<D.columns();j++) if (memo->D(i,j)!=D(i,j)) return false;
            return true;
        };
        std::vector<rvec_t> rho;
        if (sameD()) rho = memo->rho;
        else
        {
            rho = lat->CollocateDensity(D, phase, A, N_L, ecut_L);
            memo->D=D; memo->rho=rho; memo->valid=true;
        }
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
// short-circuits ContractG_ERI3), so `kernel` stays empty.
G_ERI3 GPW_Evaluator::Repulsion3CTensor() const
{
    G_ERI3 g;
    g.volume=itsGrid->Volume();
    for (const ivec3_t& dm : itsGrid->Gs()) g.columns.push_back({dm,{}});
    g.apply=MakeCollocator(/*coulomb*/true);
    return g;
}
// Overlap 3-centre tensor: the same analytic map, no kernel (the density's own rho-tilde).
G_ERI3 GPW_Evaluator::Overlap3CTensor() const
{
    G_ERI3 g;
    g.volume=itsGrid->Volume();
    for (const ivec3_t& dm : itsGrid->Gs()) g.columns.push_back({dm,{}});
    g.apply=MakeCollocator(/*coulomb*/false);
    return g;
}

// The potential->KS bridge (Band_FT_IBS::MakeOverlap / isPW_DFT_Evaluator) -- the ANALYTIC integrate-back,
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
    return itsLat->IntegratePotential(V_L, CellPhase(), A, itsLevelN, itsLevelEcut, 1.0, screenD);
}

// The REL_CUTOFF multi-grid density-grid ladder: the fine grid (L=0, reused) plus coarser grids each a factor
// 4 lower in Ecut, down to the level that still resolves the most-diffuse orbital product (exponent 2*alpha_min).
// Built once (geometry-fixed).  Because the collocation/integrate-back are ANALYTIC (only V is ever sampled), a
// diffuse pair on its matched coarse level is accurate -- the sampling multigrid's aliasing (and its depth cap)
// are gone.
void GPW_Evaluator::EnsureLevels() const
{
    if (!itsLevels.empty()) return;
    assert(itsGrid && "GPW_Evaluator: the DFT tier requires the density grid (densityEcut!=0)");
    const double efine=itsGrid->Ecut();
    const double amax=itsLat->MaxExponent(), amin=itsLat->MinExponent();
    const double ecoarse=efine*amin/amax;             // resolves the most-diffuse pair product (exponent 2*amin)
    itsLevels.push_back(itsGrid);                     // L=0: the fine grid, reused
    double e=efine;
    while (e/4.0>=ecoarse)                            // factor-4 coarsening down to the diffuse floor
    {
        e/=4.0;
        auto g=std::make_shared<const PW_Grid_Evaluator>(itsGrid->Recip(), rvec3_t(0,0,0), e);
        // RESOLUTION GUARD: keep a level only if its grid SPACING still quadratures the SHARPEST pair product
        // the REL_CUTOFF assignment can send it, p_max = 2 alpha_max ecut_l/ecut_fine: require h <= 1/sqrt(p_max)
        // (trapezoid/Poisson error of exp(-p r^2) on spacing h is ~e^{-2(pi sigma/h)^2} <= e^{-2 pi^2} ~ 3e-9 at
        // h = sigma).  This replaces a naive minimum-N floor on BOTH sides: it REJECTS a degenerate few-point
        // grid (an FFT grid that saturates at N~1-2 no longer scales like sqrt(Ecut); a sigma~2 pair on h=7 lost
        // percent-level charge -- the failed crystal gate), and it ADMITS the properly-coarse N~3-5 grids the
        // ultra-diffuse pairs want (keeping them on a fine grid made their big-box cross-cell offset sums the
        // per-iteration hotspot).  Pairs below the last level use the coarsest surviving one (finer than needed).
        const ivec3_t Ng=g->FFTGrid();
        const double  h=itsCell.GetMaximumCellEdge()/double(std::min(std::min(Ng.x,Ng.y),Ng.z));
        const double  pmax=2.0*amax*e/efine;
        if (h*h*pmax > 1.0) break;
        itsLevels.push_back(g);
    }
    for (const auto& g : itsLevels)
    {
        itsLevelN.push_back(g->FFTGrid());
        itsLevelEcut.push_back(g->Ecut());
    }
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

// The periodic 1E matrices: analytic single lattice sums (LatticeSum1E) over the SAME image set the density
// collocation uses -- the COMPLETE-Bloch scheme, self-consistent as Rcut grows.  The overlap becomes positive-
// definite once the image sphere is large enough (a truncated single sum can be indefinite; overlap integrals
// are cheap, so a generous Rcut is the fix).  KineticMatrix is <p^2> (no 1/2 -- the Hamiltonian applies it).
chmat_t GPW_Evaluator::OverlapMatrix()                 const {return itsLat->MakeOverlap(itsR,itsPhase);}
chmat_t GPW_Evaluator::KineticMatrix()                 const {return itsLat->MakeKinetic(itsR,itsPhase);}
chmat_t GPW_Evaluator::NuclearMatrix(const Structure* cl) const {return itsLat->MakeNuclear(itsR,itsPhase,cl);}

// The PP-quadrature integration mesh: a uniform lattice mesh whose Nyquist resolution follows the density
// cutoff (CreateIntegrationMesh's "GPW / Nyquist path").  The DFT tier (density grid) must be on: PP assembly
// only ever runs inside an SCF, which needs the density collocation grid anyway.
qcMesh::MeshParams GPW_Evaluator::PPMeshParams() const
{
    assert(itsGrid && "GPW_Evaluator: the external PP requires the DFT density grid (densityEcut!=0: <0 auto, >0 explicit)");
    qcMesh::MeshParams mp;
    mp.eCut=itsGrid->Ecut();
    return mp;
}

// Local PP: assembled in G-SPACE from the analytic form factor, IDENTICALLY to PW_Evaluator::LocalPotential-
// Matrix -- Vtilde(dG) = (1/Omega) Sum_a v_loc(Z_a,|dG|^2) e^{-i dG.tau_a}, dG=0 DROPPED, then <chi|Vtilde|chi>
// via OverlapMatrix (the collocation adjoint reconstructs V_loc(r) on the density grid and quadratures it).
// This inherits the PW G=0 / FormFactorG0-alignment convention EXACTLY, so the energy is box-independent (a
// real-space quadrature of the raw -Zion/r-tailed V_loc has a cell-size-dependent mean -> box drift).  The KB
// nonlocal (localized, no Coulomb tail, no G=0 issue) stays real-space (MakeSeparablePP).
chmat_t GPW_Evaluator::MakeLocalPP(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    assert(itsGrid && "GPW_Evaluator: the local PP needs the density grid (densityEcut!=0: <0 auto, >0 explicit)");
    const UnitCell& B=itsGrid->Recip().GetCell();
    // ANALYTIC integrate-back of the G-space local-PP form factor, on the ladder with a STIFFENED pair->level
    // requirement (relCutoffScale=6).  The energy error of a (pair x field) product decays with the SUM of the
    // two spectra: V_loc is spectrally BROAD (unlike the smooth V_H/V_xc), so each pair must sit ~6x finer than
    // the smooth-field calibration -- at the smooth setting the mid pairs' e^{-2.5} tails against the deep well
    // were a Ha-scale term (Si SR Gamma -10.72 vs CP2K -7.115 with charge exact).  Fine-only (every pair on the
    // fine grid) is exact but UNTENABLE at a large fine grid (NaF's auto-160: an ultra-diffuse pair's box on
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
        for (Atom* a : *cl) acc += loc.FormFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        return acc/itsGrid->Volume();
    };
    EnsureLevels();
    const size_t K=itsLevels.size();
    std::vector<rvec_t> V_L(K);
    for (size_t L=0;L<K;L++)
    {
        ΔG_Map vmapL;
        for (const ivec3_t& dm : itsLevels[L]->Gs()) vmapL[dm]=Vt(dm);   // restrict to level L's {G}
        V_L[L]=itsLevels[L]->RhoOnGrid(vmapL);
    }
    return itsLat->IntegratePotential(V_L, CellPhase(), itsCell, itsLevelN, itsLevelEcut,
                                      /*relCutoffScale*/6.0);            // sharp-field level assignment
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
        const size_t n=itsN;
        std::vector<rvec3_t> Rs(itsRc.size());
        cvec_t ph(itsRc.size());
        for (size_t r=0;r<itsRc.size();r++) { Rs[r]=-1.0*itsRc[r]; ph[r]=std::conj(itsPhaseC[r]); }
        mat_t<dcmplx> V(n,n,dcmplx(0.0));
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
                        cvec_t bt=itsLat->MakeOverlap(Rs, ph, g);
                        for (size_t i=0;i<n;i++) b[i]+=bt[i];
                    }
                    for (size_t i=0;i<n;i++)
                        for (size_t j=0;j<n;j++) V(i,j)+=D*b[i]*std::conj(b[j]);
                }
            }
        }
        chmat_t H=blazem::zeroH<dcmplx>(n);
        for (size_t i=0;i<n;i++)
        {
            H(i,i)=dcmplx(std::real(V(i,i)),0.0);
            for (size_t j=i+1;j<n;j++) H(i,j)=V(i,j);
        }
        return H;
    }

    qcMesh::Mesh mesh=cl->CreateIntegrationMesh(PPMeshParams());
    const rvec3vec_t& R=mesh.Points();
    const rvec_t&     W=mesh.Weights();
    size_t n=itsN, npts=mesh.size();

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

    mat_t<dcmplx> V(n,n,dcmplx(0.0));
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
                for (size_t i=0;i<n;i++)
                    for (size_t j=0;j<n;j++) V(i,j)+=D*b[i]*std::conj(b[j]);
            }
        }
    }
    chmat_t H=blazem::zeroH<dcmplx>(n);            // project to Hermitian (upper triangle; diagonal real)
    for (size_t i=0;i<n;i++)
    {
        H(i,i)=dcmplx(std::real(V(i,i)),0.0);
        for (size_t j=i+1;j<n;j++) H(i,j)=V(i,j);
    }
    return H;
}

// Cache key: the molecular basis's geometry-aware ID pins the radials + centres (so the cell geometry is in
// here via the atom positions); k + the translation count distinguish the periodic block.
std::string GPW_Evaluator::IDFragment() const
{
    // Include the density-grid cutoff: the collocation tensor (Repulsion3C/Overlap3C) is built on that grid, so
    // the framework cache (keyed by BasisSetID) must distinguish GPW bases that differ only in densityEcut.
    return "|mol="+itsOrb->BasisSetID()
         +"|k="+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
         +"|nR="+std::to_string(itsR.size())+"|nRc="+std::to_string(itsRc.size())   // both image sets
         +"|dEcut="+(itsGrid?std::to_string(itsGrid->Ecut()):std::string("0"));
}

} //namespace
