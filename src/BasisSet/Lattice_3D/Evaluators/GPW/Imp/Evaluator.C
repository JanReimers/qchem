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
#include <algorithm>
#include <cblas.h>   // cblas_zgemm: the OverlapMatrix integrate-back contraction (NaF profile hotspot ~44%) --
                     // OpenBLAS's SIMD/threaded complex GEMM.  blaze's own product is scalar in our TUs
                     // (BLAZE_BLAS_MODE=0 -- Blaze_Import's define rides the Blaze cmake target we don't link).
module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.Blaze;       // rvec_t, rmat_t, rsmat_t, blazem::zeroH<dcmplx>
import qchem.Vector3D;    // vec3_t + rvec3_t / rvec3vec_t arithmetic (r - R, componentwise add)
import qchem.Mesh.Quadrature; // qcMesh::WeightedOverlap / Overlap (the real-space PP quadrature primitives)
import qchem.ScalarFunction;  // ScalarFunction<double> (the V_loc / beta*Ylm fields handed to the quadrature)
import qchem.Math;            // norm, Pi, sqrt (the real spherical harmonics for the KB projectors)
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
} //anon

GPW_Evaluator::GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                             double densityEcut, const rvec3_t& kFrac, double Rcut, double collRcut,
                             double cutoffFactor)
    : itsMol(std::move(mol))
    , itsk(kFrac)
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

    // Two decoupled {R}+{e^{ik.R}} weighted point sets: the OVERLAP/1E set (Rcut -- must be large enough that
    // the analytic single-sum overlap is positive-definite) and the COLLOCATION set (collRcut -- the local
    // density's orbital reach, much smaller; the collocation re-sums images at every grid point every SCF
    // iteration, so shrinking it is the key multi-k speed-up).  CellsInSphere is inversion-symmetric, so the
    // lattice-sum matrices come out Hermitian (real at Gamma, where every phase is 1).  collRcut<=0 reuses the
    // overlap set (backward-compatible: Gamma/finite runs collocate on the same set).
    BuildImages(cell, Rcut, itsk, itsR, itsPhase);
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

// chi_i(r_g) on the density grid (Bloch-summed Gaussians; real at Gamma).  CACHED in itsPhiOnGrid: it is a
// pure function of geometry (grid points + orbitals + image set), so it is identical across SCF iterations,
// and both BuildWeights (one-time tensor) and the PER-ITERATION integrate-back OverlapMatrix reuse it.  The
// Bloch sum is O(Npts x |image set|) and DENSITY-INDEPENDENT -- recomputing it every OverlapMatrix call was the
// dominant cost at a large overlap Rcut (many images), so it is built once here.
const mat_t<dcmplx>& GPW_Evaluator::PhiOnGrid() const
{
    assert(itsGrid && "GPW_Evaluator: DFT tier requires a density grid (densityEcut!=0: <0 auto, >0 explicit)");
    if (itsPhiOnGrid.rows()!=0) return itsPhiOnGrid;   // built once, reused
    const rvec3vec_t& pts=itsGrid->GridPoints();
    size_t Npts=pts.size(), n=itsN;
    itsPhiOnGrid.resize(Npts,n);
    for (size_t g=0; g<Npts; g++)
    {
        cvec_t v=Eval(pts[g]);      // chi_i^k(r_g) = Sum_R e^{ik.R} chi_i(r_g - R) (real at Gamma)
        for (size_t i=0;i<n;i++) itsPhiOnGrid(g,i)=v[i];
    }
    return itsPhiOnGrid;
}

// The D-free 3-centre weight tensor W_c(i,j) = (1/Omega) integral chi_i chi_j e^{-iG_c.r} =
// GridCoeff(ForwardFFT(chi_i chi_j on grid), G_c) -- one FFT per orbital pair, one column per grid {G}.  NO
// kernel (the overlap-metric weight); Repulsion3CTensor adds the Poisson kernel.  The framework caches the
// tensors this feeds, so this is a plain stateless build.
G_ERI3 GPW_Evaluator::BuildWeights() const
{
    const mat_t<dcmplx>& Phi=PhiOnGrid();
    size_t Npts=Phi.rows(), n=itsN;
    const std::vector<ivec3_t>& Gs=itsGrid->Gs();
    G_ERI3 g;
    g.volume=itsGrid->Volume();
    g.columns.resize(Gs.size());
    g.weights.assign(Gs.size(), mat_t<dcmplx>(n,n,dcmplx(0.0)));
    for (size_t c=0;c<Gs.size();c++) g.columns[c].dm=Gs[c];

    // W_c(i,j) = (1/Omega) integral chi_i^k conj(chi_j^k) e^{-iG_c.r}: the density-convention weight.  The KET
    // (j) slot is CONJUGATED, matching the PHYSICAL density rho = Sum_ij D_ij chi_i conj(chi_j) (= Sum_occ
    // |psi|^2 with psi=Sum c_i chi_i, D_ij=Sum_occ c_i conj(c_j)) -- the SAME convention as IrrepCD::operator()
    // (trans(phi) D conj(phi)) and the plane-wave delta path (rho-tilde(G)=Sum_{G_i-G_j=G} D_ij).  ContractG_ERI3
    // forms Sum_ij W_c(i,j) D_ij, so this yields the physical rho-tilde.  Conjugating the BRA (i) slot instead
    // gives Sum_ij D_ij conj(chi_i) chi_j = the TRANSPOSE-density D^T -- a DIFFERENT real field at genuinely
    // complex k (it over-binds Hartree/XC; the k=1/4 complex-k fix, doc/GPWPlan.md).  At Gamma / half-
    // integer k the orbitals are real so conj is a no-op and this is byte-identical to the old form.
    // NOT (i,j)-symmetric at general k, so loop the full n^2.
    cvec_t prod(Npts);
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            for (size_t p=0; p<Npts; p++) prod[p]=Phi(p,i)*std::conj(Phi(p,j));
            cvec_t P=itsGrid->ForwardFFT(prod);          // complex-input FFT
            for (size_t c=0;c<Gs.size();c++)
                g.weights[c](i,j)=itsGrid->GridCoeff(P, Gs[c]);
        }
    return g;
}

// The matrix-free density->rho-tilde / ->V_H map (the G_ERI3::apply realization, CP2K analytic collocation).
// Collocates rho = Sum_ij D_ij chi_i chi_j analytically on compact boxes (LatticeSum1E::CollocateDensity,
// modulo-wrapped -- NO image sum, NO Rcut, NO ringing) then FFTs to G-space.  rho-tilde(G)=FFT[rho]/Npts (the
// grid quadrature's 1/Omega and Omega/Npts cancel).  For Coulomb the diagonal Poisson kernel 4pi/|G|^2 is folded
// in (V_H = 4pi rho-tilde/G^2; G=0 -> 0).  REPLACES the dense O(N_G n^2) W-tensor (BuildWeights) -- same map,
// applied matrix-free.  Gamma: the block density is real (Re(D)).  The closure keeps the FFT grid + molecular
// basis alive (captured shared_ptr) since it lives in the framework-cached G_ERI3.
std::function<ΔG_Map(const chmat_t&)> GPW_Evaluator::MakeCollocator(bool coulomb) const
{
    const UnitCell A = itsGrid->Recip().GetCell().MakeReciprocalCell();   // direct cell (by value)
    const ivec3_t  N = itsGrid->FFTGrid();
    auto grid = itsGrid;                                    // shared_ptr -> the closure keeps the grid alive
    auto mol  = itsMol;                                     // shared_ptr -> keep the molecular basis (lat) alive
    const Molecule::LatticeSum1E* lat = itsLat;
    std::vector<ivec3_t> Gs = itsGrid->Gs();
    ReciprocalLattice recip = itsGrid->Recip();
    return [A,N,grid,mol,lat,Gs,recip,coulomb](const chmat_t& D) -> ΔG_Map
    {
        const size_t n=D.rows();
        rmat_t Dr(n,n);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) Dr(i,j)=std::real(dcmplx(D(i,j)));  // Gamma: real block D
        rvec_t rho = lat->CollocateDensity(Dr, A, N);       // analytic, modulo-wrapped collocation
        cvec_t rhoTilde = grid->ForwardFFT(rho);            // rho-tilde(G)
        ΔG_Map out;
        for (const ivec3_t& dm : Gs)
        {
            const dcmplx c = grid->GridCoeff(rhoTilde, dm);
            out[dm] = coulomb ? recip.CoulombKernel(dm)*c : c;
        }
        return out;
    };
}
// Coulomb 3-centre tensor: the single-r collocation weights + the diagonal Poisson kernel 4pi/|G_c|^2.
// NOTE: still the SAMPLING (BuildWeights) path -- the analytic MakeCollocator closure is BUILT + validated
// (GPW.AnalyticCollocationCrystalChargeConservation) but NOT wired in yet: single-grid analytic collocation is
// impractically slow for a full SCF (diffuse SIPP products x the image sum), so it awaits the REL_CUTOFF
// multigrid (Increment D) that makes it fast.  Then Repulsion/Overlap3CTensor switch to g.apply=MakeCollocator.
G_ERI3 GPW_Evaluator::Repulsion3CTensor() const
{
    G_ERI3 g=BuildWeights();
    const ReciprocalLattice& recip=itsGrid->Recip();
    g.kernel.resize(g.columns.size());
    for (size_t c=0;c<g.columns.size();c++) g.kernel[c]=recip.CoulombKernel(g.columns[c].dm);
    return g;
}
// Overlap 3-centre tensor: the same single-r weights, empty kernel.
G_ERI3 GPW_Evaluator::Overlap3CTensor() const {return BuildWeights();}

// The potential->KS-matrix bridge (collocation's adjoint): inverse-FFT Vtilde over the density grid to V(r),
// then <chi_i|V|chi_j> = integral chi_i V chi_j by grid quadrature.  Vtilde is sampled over the grid's own {G}
// (which matches the fit basis GPW created, so a Hartree/XC Vtilde covers exactly these).
// The dynamic potential->KS bridge.  With the multi-grid path ON (itsUseMG), the DYNAMIC (per-iteration) smooth
// Hartree+XC integrate-back routes through the level ladder; the sharp STATIC local PP stays on the dense fine
// path (MakeLocalPP calls DenseOverlapMatrix directly) -- coarsening a diffuse pair against the sharp external
// PP is catastrophic (Si Gamma -21.4), but the smooth V_H+V_xc coarsen cleanly (Si Gamma -8.2485 vs -8.24758,
// grid tolerance).  Default OFF (dense) so committed anchors are byte-identical until the win is NaF-validated.
// The potential->KS bridge (Band_FT_IBS::MakeOverlap / isPW_DFT_Evaluator).  STILL the sampling path (dispatch
// dense/multigrid): the ANALYTIC integrate-back (LatticeSum1E::IntegratePotential) is built + validated but not
// wired in -- it awaits the multigrid (Increment D), same as the density side (see Repulsion3CTensor).
chmat_t GPW_Evaluator::OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    if (itsUseMG) return MultiGridOverlapMatrix(Vtilde);
    return DenseOverlapMatrix(Vtilde);
}
chmat_t GPW_Evaluator::DenseOverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    const mat_t<dcmplx>& Phi=PhiOnGrid();
    ΔG_Map vmap;
    for (const ivec3_t& dm : itsGrid->Gs()) vmap[dm]=Vtilde(dm);
    rvec_t V=itsGrid->RhoOnGrid(vmap);          // V(r) on the grid (real: Hartree/XC/local-PP are real fields)
    size_t Npts=Phi.rows(), n=itsN;
    const double w=itsGrid->Volume()/double(Npts);   // uniform quadrature weight (== Integral's Omega/Npts)
    // <chi_i^k|V|chi_j^k> = w Sum_p conj(Phi(p,i)) V[p] Phi(p,j) = w (Phi^H (V .* Phi))(i,j).  Form (V .* Phi)
    // once (one streaming pass over Phi), then Phi^H (V .* Phi) via OpenBLAS zgemm -- a genuinely SIMD/threaded
    // complex GEMM.  This was THE per-iteration hotspot in the NaF profile (~44%): the plain triple loop re-read
    // the ~130 MB Phi (Npts x n) n^2/2 times (memory-bound) AND ran scalar; the GEMM streams Phi ~once and
    // vectorizes.  (blaze's own product gave NO speedup here -- BLAZE_BLAS_MODE=0 in our TUs, so it too is scalar
    // complex; hence the direct cblas call.)  V is real -> M Hermitian; project to the upper triangle (real
    // diagonal).  Translation invariance is carried by the Bloch orbital in PhiOnGrid (Rcut>=2a; see the plan).
    itsVPhiBuf.resize(Npts,n);                        // reused scratch: allocates only on the first call (dims are
    mat_t<dcmplx>& VPhi=itsVPhiBuf;                   // geometry-fixed), no-op resize thereafter -- avoids ~130 MB
    for (size_t j=0;j<n;j++)                          // alloc/free per call.  Column-major Phi: p (row) inner is stride-1.
        for (size_t p=0;p<Npts;p++) VPhi(p,j)=V[p]*Phi(p,j);
    mat_t<dcmplx> C(n,n,dcmplx(0.0));                 // C = w Phi^H VPhi  (n x n; tiny, leave as a local)
    const dcmplx alpha(w,0.0), beta(0.0,0.0);
    cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
                (blasint)n, (blasint)n, (blasint)Npts,
                &alpha, Phi.data(),  (blasint)Phi.spacing(),
                        VPhi.data(), (blasint)VPhi.spacing(),
                &beta,  C.data(),    (blasint)C.spacing());
    chmat_t M=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
    {
        M(i,i)=dcmplx(std::real(C(i,i)),0.0);
        for (size_t j=i+1;j<n;j++) M(i,j)=C(i,j);
    }
    return M;
}

// The PATCHED integrate-back: build V(r) exactly as above (inverse-FFT of Vtilde over the density grid), then
// delegate the <chi_i|V|chi_j> quadrature to the molecular side, which contracts each pair only on the overlap
// of its two orbitals' Gaussian supports (the sub-eps grid points are dropped -> bit-consistent with the dense
// GEMM above to the screening tolerance).  Uses the COLLOCATION image set (itsRc/itsPhaseC) -- the same Bloch
// sum PhiOnGrid/Eval use -- so the two paths sample the identical chi_i^k(r_g).  Opt-in scaffold; see the header.
chmat_t GPW_Evaluator::PatchedOverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    ΔG_Map vmap;
    for (const ivec3_t& dm : itsGrid->Gs()) vmap[dm]=Vtilde(dm);
    rvec_t V=itsGrid->RhoOnGrid(vmap);                       // V(r) on the density grid (real)
    const rvec3vec_t& pts=itsGrid->GridPoints();
    const double w=itsGrid->Volume()/double(V.size());       // uniform quadrature weight Omega/Npts
    return itsLat->MakePotentialMatrix(pts, itsRc, itsPhaseC, V, w);
}

// The multi-grid density-grid ladder: the fine grid (L=0, reused) plus coarser grids each a factor 4 lower in
// Ecut, down to the level that still resolves the most-diffuse orbital product (exponent 2*alpha_min).  Built
// once (geometry-fixed) and cached with each level's points / cutoff / quadrature weight.
void GPW_Evaluator::EnsureLevels() const
{
    if (!itsLevels.empty()) return;
    assert(itsGrid && "GPW_Evaluator: the multi-grid integrate-back requires the density grid (densityEcut!=0)");
    const double efine=itsGrid->Ecut();
    const double amax=itsLat->MaxExponent(), amin=itsLat->MinExponent();
    const double ecoarse=efine*amin/amax;             // resolves the most-diffuse pair product (exponent 2*amin)
    itsLevels.push_back(itsGrid);                     // L=0: the fine grid, reused
    double e=efine;
    const size_t maxL=itsMGMaxLevels>0 ? size_t(itsMGMaxLevels) : 8;  // depth cap (REL_CUTOFF-style safety)
    while (e/4.0>=ecoarse && itsLevels.size()<maxL)   // factor-4 coarsening; cap the ladder depth
    {
        e/=4.0;
        itsLevels.push_back(std::make_shared<const PW_Grid_Evaluator>(itsGrid->Recip(), rvec3_t(0,0,0), e));
    }
    for (const auto& g : itsLevels)
    {
        itsLevelPts.push_back(g->GridPoints());
        itsLevelEcut.push_back(g->Ecut());
        itsLevelW.push_back(g->Volume()/double(g->GridPoints().size()));
    }
}

// The MULTI-GRID integrate-back: restrict Vtilde to each level's own {G} (the low-pass) -> V(r) on that level,
// then let the molecular side contract each orbital pair on the coarsest level resolving its product.  The
// per-(level,orbital) chi columns are cached molecular-side, so per-iteration only V_L is rebuilt here.
chmat_t GPW_Evaluator::MultiGridOverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    EnsureLevels();
    const size_t K=itsLevels.size();
    std::vector<rvec_t> V_L(K);
    for (size_t L=0;L<K;L++)
    {
        ΔG_Map vmapL;
        for (const ivec3_t& dm : itsLevels[L]->Gs()) vmapL[dm]=Vtilde(dm);   // restrict to level L's {G} (low-pass)
        V_L[L]=itsLevels[L]->RhoOnGrid(vmapL);
    }
    return itsLat->MakePotentialMatrixMG(itsLevelPts, itsLevelEcut, itsRc, itsPhaseC, V_L, itsLevelW);
}

// Bloch sum of the Gaussian orbitals, chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R), over the COLLOCATION set (the
// orbital reach; == the overlap set unless collRcut decoupled it).  At Gamma every phase is 1, so the
// imaginary part is exactly zero (the sum reduces to the real molecular sum).
cvec_t GPW_Evaluator::Eval(const rvec3_t& r) const
{
    cvec_t v(itsN, dcmplx(0.0));
    for (size_t k=0;k<itsRc.size();k++)
    {
        rvec_t chi=(*itsOrb)(r-itsRc[k]);
        for (size_t i=0;i<itsN;i++) v[i]+=itsPhaseC[k]*chi[i];
    }
    return v;
}

cvec3vec_t GPW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    cvec3vec_t v(itsN, vec3_t<dcmplx>(dcmplx(0.0),dcmplx(0.0),dcmplx(0.0)));
    for (size_t k=0;k<itsRc.size();k++)
    {
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
    // DENSE (sampling) integrate-back of the G-space local-PP form factor (still the committed path; the analytic
    // OverlapMatrix awaits the multigrid, Increment D).
    return DenseOverlapMatrix([&](const ivec3_t& dm)->dcmplx
    {
        if (dm.x==0 && dm.y==0 && dm.z==0) return dcmplx(0.0);        // drop dG=0 (alignment carries it)
        rvec3_t dG=B.ToCartesian(rvec3_t(dm));
        double  g2=dG*dG;
        dcmplx  acc(0.0);                                            // form factor x structure factor
        for (Atom* a : *cl) acc += loc.FormFactor(a->itsZ,g2)*std::exp(dcmplx(0.0,-(dG*a->itsR)));
        return acc/itsGrid->Volume();
    });
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
