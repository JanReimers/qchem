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

namespace qchem::BasisSet::Lattice_3D
{

namespace
{
// Widen a real symmetric lattice-sum matrix to the complex Hermitian type the concept/mixin speak.  At the
// Gamma point the lattice sums are real, so the imaginary part is exactly zero (asserted by the tests).
chmat_t Complexify(const rsmat_t& S)
{
    size_t n=S.rows();
    chmat_t C=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            C(i,j)=dcmplx(S(i,j),0.0);
    return C;
}
} //anon

GPW_Evaluator::GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                             double densityEcut, const rvec3_t& kFrac, double Rcut)
    : itsMol(std::move(mol))
    , itsk(kFrac)
{
    assert(std::fabs(itsk.x)+std::fabs(itsk.y)+std::fabs(itsk.z) < 1e-12
           && "GPW_Evaluator: only the Gamma point (k=0) is supported this increment");

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

    // The real-space translation set {R} (Cartesian).  Rcut<=0 keeps only the home cell (R=0 == the finite
    // molecule); a positive Rcut adds the images inside the sphere.  CellsInSphere yields an inversion-
    // symmetric set, so the lattice-sum matrices come out symmetric.
    if (Rcut>0.0)
        for (const auto& n : cell.CellsInSphere(Rcut))
            itsR.push_back(cell.ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z))));
    else
        itsR.push_back(rvec3_t(0,0,0));
    assert(!itsR.empty());

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
rmat_t GPW_Evaluator::PhiOnGrid() const
{
    assert(itsGrid && "GPW_Evaluator: DFT tier requires a density grid (construct with densityEcut>0)");
    const rvec3vec_t& pts=itsGrid->GridPoints();
    size_t Npts=pts.size(), n=itsN;
    rmat_t Phi(Npts,n);
    for (size_t g=0; g<Npts; g++)
    {
        cvec_t v=Eval(pts[g]);
        for (size_t i=0;i<n;i++) Phi(g,i)=std::real(v[i]);
    }
    return Phi;
}

// The D-free 3-centre weight tensor W_c(i,j) = (1/Omega) integral chi_i chi_j e^{-iG_c.r} =
// GridCoeff(ForwardFFT(chi_i chi_j on grid), G_c) -- one FFT per orbital pair, one column per grid {G}.  NO
// kernel (the overlap-metric weight); Repulsion3CTensor adds the Poisson kernel.  The framework caches the
// tensors this feeds, so this is a plain stateless build.
G_ERI3 GPW_Evaluator::BuildWeights() const
{
    rmat_t Phi=PhiOnGrid();
    size_t Npts=Phi.rows(), n=itsN;
    const std::vector<ivec3_t>& Gs=itsGrid->Gs();
    G_ERI3 g;
    g.volume=itsGrid->Volume();
    g.columns.resize(Gs.size());
    g.weights.assign(Gs.size(), mat_t<dcmplx>(n,n,dcmplx(0.0)));
    for (size_t c=0;c<Gs.size();c++) g.columns[c].dm=Gs[c];

    rvec_t prod(Npts);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            for (size_t p=0; p<Npts; p++) prod[p]=Phi(p,i)*Phi(p,j);   // chi_i chi_j on the grid (real at Gamma)
            cvec_t P=itsGrid->ForwardFFT(prod);
            for (size_t c=0;c<Gs.size();c++)
            {
                dcmplx w=itsGrid->GridCoeff(P, Gs[c]);   // (1/Omega) integral chi_i chi_j e^{-iG_c.r} (ONE r)
                g.weights[c](i,j)=w;
                g.weights[c](j,i)=w;                     // symmetric in (i,j)
            }
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
    rmat_t Phi=PhiOnGrid();
    ΔG_Map vmap;
    for (const ivec3_t& dm : itsGrid->Gs()) vmap[dm]=Vtilde(dm);
    rvec_t V=itsGrid->RhoOnGrid(vmap);          // V(r) on the grid (real)
    size_t Npts=Phi.rows(), n=itsN;
    rvec_t f(Npts);
    chmat_t M=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            for (size_t p=0; p<Npts; p++) f[p]=Phi(p,i)*V[p]*Phi(p,j);
            M(i,j)=dcmplx(itsGrid->Integral(f), 0.0);   // real at Gamma
        }
    return M;
}

// Bloch sum of the Gaussian orbitals, chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R).  At Gamma every phase is 1, so
// this is the real sum of the molecular basis evaluated at the image-shifted points (widened to complex).
cvec_t GPW_Evaluator::Eval(const rvec3_t& r) const
{
    rvec_t acc=(*itsOrb)(r-itsR[0]);
    for (size_t k=1;k<itsR.size();k++) acc += (*itsOrb)(r-itsR[k]);
    cvec_t v(itsN);
    for (size_t i=0;i<itsN;i++) v[i]=dcmplx(acc[i],0.0);
    return v;
}

cvec3vec_t GPW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    rvec3vec_t acc=itsOrb->Gradient(r-itsR[0]);
    for (size_t k=1;k<itsR.size();k++)
    {
        rvec3vec_t g=itsOrb->Gradient(r-itsR[k]);
        for (size_t i=0;i<itsN;i++) acc[i]=acc[i]+g[i];
    }
    cvec3vec_t v(itsN);
    for (size_t i=0;i<itsN;i++) v[i]=vec3_t<dcmplx>(dcmplx(acc[i].x,0.0),dcmplx(acc[i].y,0.0),dcmplx(acc[i].z,0.0));
    return v;
}

// The periodic 1E matrices: delegate the lattice sum to the molecular basis (it owns the Gaussian kernels),
// then widen to complex.  KineticMatrix is <p^2> (no 1/2 -- matches PW_Evaluator; the Hamiltonian applies it).
chmat_t GPW_Evaluator::OverlapMatrix()                 const {return Complexify(itsLat->MakeOverlap(itsR));}
chmat_t GPW_Evaluator::KineticMatrix()                 const {return Complexify(itsLat->MakeKinetic(itsR));}
chmat_t GPW_Evaluator::NuclearMatrix(const Structure* cl) const {return Complexify(itsLat->MakeNuclear(itsR,cl));}

// Cache key: the molecular basis's geometry-aware ID pins the radials + centres (so the cell geometry is in
// here via the atom positions); k + the translation count distinguish the periodic block.
std::string GPW_Evaluator::IDFragment() const
{
    // Include the density-grid cutoff: the collocation tensor (Repulsion3C/Overlap3C) is built on that grid, so
    // the framework cache (keyed by BasisSetID) must distinguish GPW bases that differ only in densityEcut.
    return "|mol="+itsOrb->BasisSetID()
         +"|k="+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
         +"|nR="+std::to_string(itsR.size())
         +"|dEcut="+(itsGrid?std::to_string(itsGrid->Ecut()):std::string("0"));
}

} //namespace
