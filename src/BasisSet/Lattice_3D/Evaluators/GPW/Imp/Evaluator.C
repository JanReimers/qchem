// File: BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C  GPW_Evaluator implementation.
module;
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.Blaze;       // rvec_t, rmat_t, rsmat_t, blazem::zeroH<dcmplx>
import qchem.Vector3D;    // vec3_t + rvec3_t / rvec3vec_t arithmetic (r - R, componentwise add)
import qchem.Mesh.Quadrature; // qcMesh::WeightedOverlap / Overlap (the real-space PP quadrature primitives)
import qchem.ScalarFunction;  // ScalarFunction<double> (the V_loc / beta*Ylm fields handed to the quadrature)
import qchem.Math;            // norm, Pi, sqrt (the real spherical harmonics for the KB projectors)
import qchem.Structure;       // Structure / Atom (the PP centres + Z, and CreateIntegrationMesh)

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
} //anon

GPW_Evaluator::GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                             double densityEcut, const rvec3_t& kFrac, double Rcut, double collRcut)
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

    // The DFT tier's density/collocation grid: a plane-wave grid at Gamma (the density is cell-periodic
    // whatever the orbital k) and the caller's density cutoff.  The fit basis GPW_IBS creates is built over
    // THIS grid, so rho-tilde's {G} matches the fitter's.  Off (null) when no DFT is needed (1E-only tests).
    if (densityEcut>0.0)
        itsGrid=std::make_shared<const PW_Grid_Evaluator>(
                    ReciprocalLattice(cell.MakeReciprocalCell()), rvec3_t(0,0,0), densityEcut);
}

// chi_i(r_g) on the density grid (Bloch-summed Gaussians; real at Gamma).  Computed on demand -- NOT cached in
// a member: the framework caches the DERIVED tensor (Repulsion3C/Overlap3C via Band_FT_IBS's DB_Cache), so this
// runs a small, bounded number of times (once per one-time tensor build; and per OverlapMatrix call, whose
// per-iteration quadrature dominates it anyway).  Keeping the evaluator stateless keeps the caching in one place.
mat_t<dcmplx> GPW_Evaluator::PhiOnGrid() const
{
    assert(itsGrid && "GPW_Evaluator: DFT tier requires a density grid (construct with densityEcut>0)");
    const rvec3vec_t& pts=itsGrid->GridPoints();
    size_t Npts=pts.size(), n=itsN;
    mat_t<dcmplx> Phi(Npts,n);
    for (size_t g=0; g<Npts; g++)
    {
        cvec_t v=Eval(pts[g]);      // chi_i^k(r_g) = Sum_R e^{ik.R} chi_i(r_g - R) (real at Gamma)
        for (size_t i=0;i<n;i++) Phi(g,i)=v[i];
    }
    return Phi;
}

// The D-free 3-centre weight tensor W_c(i,j) = (1/Omega) integral chi_i chi_j e^{-iG_c.r} =
// GridCoeff(ForwardFFT(chi_i chi_j on grid), G_c) -- one FFT per orbital pair, one column per grid {G}.  NO
// kernel (the overlap-metric weight); Repulsion3CTensor adds the Poisson kernel.  The framework caches the
// tensors this feeds, so this is a plain stateless build.
G_ERI3 GPW_Evaluator::BuildWeights() const
{
    mat_t<dcmplx> Phi=PhiOnGrid();
    size_t Npts=Phi.rows(), n=itsN;
    const std::vector<ivec3_t>& Gs=itsGrid->Gs();
    G_ERI3 g;
    g.volume=itsGrid->Volume();
    g.columns.resize(Gs.size());
    g.weights.assign(Gs.size(), mat_t<dcmplx>(n,n,dcmplx(0.0)));
    for (size_t c=0;c<Gs.size();c++) g.columns[c].dm=Gs[c];

    // W_c(i,j) = (1/Omega) integral conj(chi_i^k) chi_j^k e^{-iG_c.r}: the density-convention weight (the i slot
    // is CONJUGATED, since rho = Sum_ij D_ij conj(chi_i) chi_j and ContractG_ERI3 forms Sum_ij W_c(i,j) D_ij).
    // NOT (i,j)-symmetric at general k, so loop the full n^2 (at Gamma both real -> collapses to the old form).
    cvec_t prod(Npts);
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            for (size_t p=0; p<Npts; p++) prod[p]=std::conj(Phi(p,i))*Phi(p,j);
            cvec_t P=itsGrid->ForwardFFT(prod);          // complex-input FFT
            for (size_t c=0;c<Gs.size();c++)
                g.weights[c](i,j)=itsGrid->GridCoeff(P, Gs[c]);
        }
    return g;
}

// Coulomb 3-centre tensor: the single-r collocation weights + the diagonal Poisson kernel 4pi/|G_c|^2 (the r_2
// / 1-over-r_12 factor; G_c=0 -> 0).  The framework's Repulsion3C(c) caches this build.
G_ERI3 GPW_Evaluator::Repulsion3CTensor() const
{
    G_ERI3 g=BuildWeights();
    const ReciprocalLattice& recip=itsGrid->Recip();
    g.kernel.resize(g.columns.size());
    for (size_t c=0;c<g.columns.size();c++) g.kernel[c]=recip.CoulombKernel(g.columns[c].dm);
    return g;
}

// Overlap 3-centre tensor: the same single-r weights, empty kernel (the density's rho-tilde, no Poisson).  The
// framework's Overlap3C(c) caches this build.
G_ERI3 GPW_Evaluator::Overlap3CTensor() const {return BuildWeights();}

// The potential->KS-matrix bridge (collocation's adjoint): inverse-FFT Vtilde over the density grid to V(r),
// then <chi_i|V|chi_j> = integral chi_i V chi_j by grid quadrature.  Vtilde is sampled over the grid's own {G}
// (which matches the fit basis GPW created, so a Hartree/XC Vtilde covers exactly these).
chmat_t GPW_Evaluator::OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    mat_t<dcmplx> Phi=PhiOnGrid();
    ΔG_Map vmap;
    for (const ivec3_t& dm : itsGrid->Gs()) vmap[dm]=Vtilde(dm);
    rvec_t V=itsGrid->RhoOnGrid(vmap);          // V(r) on the grid (real: Hartree/XC/local-PP are real fields)
    size_t Npts=Phi.rows(), n=itsN;
    const double w=itsGrid->Volume()/double(Npts);   // uniform quadrature weight (== Integral's Omega/Npts)
    // <chi_i^k | V | chi_j^k> = integral conj(chi_i^k) V chi_j^k: Hermitian (V real), fill upper triangle.
    chmat_t M=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            dcmplx s(0.0);
            for (size_t p=0; p<Npts; p++) s += std::conj(Phi(p,i))*V[p]*Phi(p,j);
            s *= w;
            M(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) auto-conj
        }
    return M;
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
    assert(itsGrid && "GPW_Evaluator: the external PP requires the DFT density grid (construct with densityEcut>0)");
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
    assert(itsGrid && "GPW_Evaluator: the local PP needs the density grid (construct with densityEcut>0)");
    const UnitCell& B=itsGrid->Recip().GetCell();
    return OverlapMatrix([&](const ivec3_t& dm)->dcmplx
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
// projection vector b_i^k = <chi_i^k | beta_p Y_lm> = Sum_R e^{ik.R} <chi_i | beta at (R_a - R)> (mesh
// quadrature of the home orbital against the image-shifted projector -- reuses the REAL qcMesh::Overlap per
// image, weighted by the phase).  At Gamma with Rcut=0 (itsR={0}, phase 1) this is exactly the finite molecular
// PP_NonLocal recipe <chi_i | beta at R_a>.  V_ij = Sum D b_i conj(b_j) is Hermitian by construction.
chmat_t GPW_Evaluator::MakeSeparablePP(const Structure* cl, const Pseudopotential::SeparablePotential_R& sep) const
{
    qcMesh::Mesh mesh=cl->CreateIntegrationMesh(PPMeshParams());
    size_t n=itsN;
    mat_t<dcmplx> V(n,n,dcmplx(0.0));
    for (size_t a=0; a<cl->GetNumAtoms(); a++)
    {
        const Atom* at=(*cl)[a];
        int Z=at->itsZ;
        for (size_t p=0; p<sep.NumProjectors(Z); p++)
        {
            int    l=sep.AngularMomentum(Z,p);
            double D=sep.Coefficient    (Z,p);
            for (int m=-l; m<=l; m++)
            {
                cvec_t b(n, dcmplx(0.0));
                for (size_t r=0; r<itsRc.size(); r++)   // orbital reach: the collocation image set
                {
                    rvec_t br=qcMesh::Overlap(mesh, *itsOrb, BetaYlmField(at->itsR-itsRc[r],sep,Z,p,l,m));
                    for (size_t i=0;i<n;i++) b[i]+=itsPhaseC[r]*br[i];
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
