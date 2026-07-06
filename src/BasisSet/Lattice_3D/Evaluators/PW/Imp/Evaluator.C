// File: BasisSet/Lattice_3D/Evaluators/PW/Imp/Evaluator.C  PW_Evaluator implementation.
module;
#include <cassert>
#include <algorithm>
#include <complex>
#include <functional>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.Evaluators.PW;
import qchem.Math;               // Pi, sqrt, cos, sin, Cube
import qchem.BasisSet.Lattice_3D.Internal.GVectors;   // BuildGs
import qchem.FFT;                // NextPow2 (the XC grid geometry)
import qchem.Blaze;              // blazem::zeroH / blazem::sum
import qchem.Vector3D;           // dot product (operator*) + vector arithmetic

namespace qchem::BasisSet::Lattice_3D
{

PW_Evaluator::PW_Evaluator(const ReciprocalLattice& recip, const rvec3_t& k, double Ecut)
    : itsRecip(recip)
    , itsk(k)
    , itsEcut(Ecut)
    // |det B| = (2 pi)^3 / |det A|, so the direct cell volume V = (2 pi)^3 / V_recip.
    , itsVolume(Cube(2*Pi)/recip.GetCell().GetCellVolume())
{
    itsG = Internal::BuildGs(itsRecip, itsk, Ecut);   // { G : 1/2|k+G|^2 < Ecut }
}

rvec3_t PW_Evaluator::GetGCartesian(const ivec3_t& m) const
{
    return itsRecip.GetCell().ToCartesian(rvec3_t(m)); // B m
}

// Uniform N1xN2xN3 grid of FRACTIONAL coordinates r=(i1/N1,i2/N2,i3/N3) -- the XC real-space grid
// (uniform weight Omega/prod(N); dG.r = 2 pi dm.r_frac makes the rho(G)<->rho(r) transform a plain DFT).
std::vector<rvec3_t> PW_Evaluator::UniformGrid(const ivec3_t& n) const
{
    std::vector<rvec3_t> g;
    g.reserve(static_cast<size_t>(n.x)*n.y*n.z);
    for (int i1=0;i1<n.x;i1++)
        for (int i2=0;i2<n.y;i2++)
            for (int i3=0;i3<n.z;i3++)
                g.push_back(rvec3_t(i1/double(n.x), i2/double(n.y), i3/double(n.z)));
    return g;
}

// Grid divisions resolving the difference set without aliasing: N > 2*(2*maxComp).  Cached (depends on itsG).
ivec3_t PW_Evaluator::AutoGrid() const
{
    if (itsAutoGrid.x==0)   // not yet computed
    {
        int m=0;
        for (const ivec3_t& g : itsG)
        {
            int ax=g.x<0?-g.x:g.x, ay=g.y<0?-g.y:g.y, az=g.z<0?-g.z:g.z;
            m=std::max(m, std::max(ax, std::max(ay,az)));
        }
        int nn=4*m+1;
        itsAutoGrid=ivec3_t(nn,nn,nn);
    }
    return itsAutoGrid;
}

// AutoGrid divisions rounded up to powers of two -- radix-2 FFT for the XC route.  A larger grid still
// resolves the difference set without aliasing, so it is at worst slightly more accurate.  Cached.
ivec3_t PW_Evaluator::FFTGrid() const
{
    if (itsFFTGrid.x==0)   // not yet computed
    {
        ivec3_t a=AutoGrid();
        itsFFTGrid=ivec3_t(int(qchem::FFT::NextPow2(a.x)), int(qchem::FFT::NextPow2(a.y)), int(qchem::FFT::NextPow2(a.z)));
    }
    return itsFFTGrid;
}

// Cartesian points of the FFT grid, raster order matching RhoOnGrid/ForwardFFT: r = A (i/N), A = direct cell.
// Cached: the N-point mesh + the 3x3 cell inversion are built once, not per DoFit.
const rvec3vec_t& PW_Evaluator::GridPoints() const
{
    if (itsGridPoints.size()==0)
    {
        std::vector<rvec3_t> frac=UniformGrid(FFTGrid());
        UnitCell A=itsRecip.GetCell().MakeReciprocalCell();   // reciprocal of B == the direct cell
        itsGridPoints.resize(frac.size());
        for (size_t q=0;q<frac.size();q++) itsGridPoints[q]=A.ToCartesian(frac[q]);
    }
    return itsGridPoints;
}

// rho(r) on the FFT grid = inverse FFT of rho-tilde: rho(r_j) = Sum_dm rho-tilde(dm) e^{+i2pi dm.j/N}.
// rho-tilde is the physical coefficient (already /Omega), so the inverse FFT takes NO 1/N normalization.
rvec_t PW_Evaluator::RhoOnGrid(const ΔG_Map& rho) const
{
    ivec3_t N=FFTGrid();
    size_t Npts=size_t(N.x)*N.y*N.z;
    cvec_t g(Npts, dcmplx(0.0));
    for (const auto& kv : rho)
    {
        int i0=((kv.first.x%N.x)+N.x)%N.x, i1=((kv.first.y%N.y)+N.y)%N.y, i2=((kv.first.z%N.z)+N.z)%N.z;
        g[(size_t(i0)*N.y+i1)*N.z+i2]=kv.second;
    }
    cvec_t rr=qchem::FFT::FFT3D(g, N, +1);
    rvec_t out(Npts);
    for (size_t i=0;i<Npts;i++) out[i]=std::real(dcmplx(rr[i]));
    return out;
}

// Forward-FFT a real-space grid field to the FULL normalised G-space grid Vtilde = FFT[V]/Npts (raster order).
// Unlike a difference-filtered map, this keeps EVERY resolved frequency, so any consumer's difference set hits.
cvec_t PW_Evaluator::ForwardFFT(const rvec_t& V) const
{
    ivec3_t N=FFTGrid();
    size_t Npts=size_t(N.x)*N.y*N.z;
    assert(V.size()==Npts);
    cvec_t g(Npts, dcmplx(0.0));
    for (size_t i=0;i<Npts;i++) g[i]=dcmplx(V[i]);
    cvec_t Vt=qchem::FFT::FFT3D(g, N, -1);
    for (size_t i=0;i<Npts;i++) Vt[i]/=double(Npts);
    return Vt;
}

// Vtilde(dm) from a ForwardFFT grid: wrap dm into [0,N) per axis (negative freq -> N-|.|) and index raster.
dcmplx PW_Evaluator::GridCoeff(const cvec_t& Vt, const ivec3_t& dm) const
{
    ivec3_t N=FFTGrid();
    // Alias-free ONLY while |dm| < N/2 per axis: outside that band the wrap folds a high frequency onto dm
    // and silently returns the wrong coefficient.  Holds by construction for relCutoff>=1 (denser fit grid
    // is protective); the real future exposure is a k!=0 XC block whose orbital difference set exceeds the
    // Gamma fit grid -- this assert makes that path fail loudly instead of aliasing.
    assert(2*(dm.x<0?-dm.x:dm.x)<N.x && 2*(dm.y<0?-dm.y:dm.y)<N.y && 2*(dm.z<0?-dm.z:dm.z)<N.z
           && "GridCoeff: |dm| outside the alias-free band -- densify the fit grid (raise relCutoff)");
    int i0=((dm.x%N.x)+N.x)%N.x, i1=((dm.y%N.y)+N.y)%N.y, i2=((dm.z%N.z)+N.z)%N.z;
    return Vt[(size_t(i0)*N.y+i1)*N.z+i2];
}

// The fitted field's coefficients over THIS basis's own {G} (a GridCoeff gather) -- for op(r) evaluation.
ΔG_Map PW_Evaluator::FieldCoeffs(const cvec_t& Vt) const
{
    ΔG_Map m;
    for (const ivec3_t& G : itsG) m[G]=GridCoeff(Vt, G);
    return m;
}

// integral f d3r on the FFT grid: uniform quadrature, weight Omega/Npts.
double PW_Evaluator::Integral(const rvec_t& f) const
{
    return blazem::sum(f)*Volume()/double(f.size());
}

// G_FieldEvaluator: inverse-transform the Hermitian coefficients c(dm) to the real field
// f(r) = Re Sum_dm c(dm) e^{i(B.dm).r}, using our OWN grid engine's B (GetGCartesian).  A point evaluation
// (loops the sparse map), NOT an FFT -- the GUI / fit-residual diagnostic path.
double PW_Evaluator::EvalField(const ΔG_Map& c, const rvec3_t& r) const
{
    dcmplx s(0.0);
    for (const auto& kv : c)
    {
        rvec3_t G = GetGCartesian(kv.first);
        double  ph = G*r;
        s += kv.second * dcmplx(cos(ph), sin(ph));
    }
    return s.real();
}

// grad f(r) = Sum_dm (B.dm) (-Im[c(dm) e^{i(B.dm).r}]) -- the gradient of the real inverse transform.
rvec3_t PW_Evaluator::EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const
{
    rvec3_t g(0.0,0.0,0.0);
    for (const auto& kv : c)
    {
        rvec3_t G = GetGCartesian(kv.first);
        double  ph = G*r;
        dcmplx  ce = kv.second * dcmplx(cos(ph), sin(ph));
        g += (-ce.imag()) * G;
    }
    return g;
}

// Plane waves are orthonormal over the cell: <G|G'> = delta_{GG'}.
chmat_t PW_Evaluator::OverlapMatrix() const
{
    size_t n=size();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // hmat_t(n) does NOT zero the off-diagonals
    for (size_t i=0; i<n; i++) S(i,i)=1.0;
    return S;
}

// <p^2> = <-nabla^2> building block (NO 1/2 -- the Hamiltonian applies it).  For a plane wave
// -nabla^2 e^{i(k+G).r} = |k+G|^2 e^{i(k+G).r}, so the matrix is diagonal in |k+G|^2.
chmat_t PW_Evaluator::KineticMatrix() const
{
    const UnitCell& B=itsRecip.GetCell();
    size_t n=size();
    chmat_t S=blazem::zeroH<dcmplx>(n);   // off-diagonals are exactly zero
    for (size_t i=0; i<n; i++)
    {
        double kG=B.GetDistance(itsk+itsG[i]); // |k+G|
        S(i,i)=kG*kG;
    }
    return S;
}

// <G|V|G'> = Vtilde(m(G) - m(G')).  Fill the upper triangle; HermitianMatrix mirrors the conjugate.
chmat_t PW_Evaluator::MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const
{
    size_t n=size();
    chmat_t V=blazem::zeroH<dcmplx>(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            V(i,j)=Vtilde(itsG[i]-itsG[j]);
    return V;
}

cvec_t PW_Evaluator::Eval(const rvec3_t& r) const
{
    size_t n=size();
    double invSqrtV=1.0/sqrt(itsVolume);
    cvec_t v(n);
    for (size_t i=0; i<n; i++)
    {
        double phase=GetGCartesian(itsG[i])*r + itsRecip.GetCell().ToCartesian(itsk)*r; // (k+G).r
        v[i]=dcmplx(cos(phase),sin(phase))*invSqrtV;
    }
    return v;
}

cvec3vec_t PW_Evaluator::EvalGradient(const rvec3_t& r) const
{
    // grad e^{i(k+G).r}/sqrt(V) = i(k+G) e^{i(k+G).r}/sqrt(V).
    const dcmplx im(0.0,1.0);
    size_t n=size();
    double invSqrtV=1.0/sqrt(itsVolume);
    rvec3_t kCart=itsRecip.GetCell().ToCartesian(itsk);
    cvec3vec_t g(n);
    for (size_t i=0; i<n; i++)
    {
        rvec3_t kG=kCart+GetGCartesian(itsG[i]); // k+G (Cartesian)
        double phase=kG*r;
        dcmplx val=dcmplx(cos(phase),sin(phase))*invSqrtV;
        g[i]=vec3_t<dcmplx>(im*kG.x*val, im*kG.y*val, im*kG.z*val);
    }
    return g;
}

std::string PW_Evaluator::IDFragment() const
{
    return std::string("|k=")+std::to_string(itsk.x)+","+std::to_string(itsk.y)+","+std::to_string(itsk.z)
                 +"|Ecut="+std::to_string(itsEcut)
                 +"|nG="+std::to_string(itsG.size());
}

} //namespace
