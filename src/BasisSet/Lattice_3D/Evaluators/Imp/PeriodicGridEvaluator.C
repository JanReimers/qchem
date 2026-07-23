// File: BasisSet/Lattice_3D/Evaluators/Imp/PeriodicGridEvaluator.C  PeriodicGridEvaluator implementation.
//
// The FFT/Poisson grid code, moved verbatim off PW_Evaluator: the {r}<->{G} transforms and the real-space
// quadrature depend ONLY on (B, Omega, N), not on the orbitals -- so they live here, k-independent, and
// PW_Evaluator (and later GPW_Evaluator) forwards its grid virtuals to a held instance.
module;
#include <cassert>
#include <complex>
#include <vector>

module qchem.BasisSet.Lattice_3D.Evaluators.PeriodicGridEvaluator;
import qchem.Math;        // cos, sin (EvalField's point evaluation)
import qchem.FFT;         // FFT3D (RhoOnGrid / ForwardFFT)
import qchem.Blaze;       // blazem::sum (Integral)
import qchem.Vector3D;    // dot product (operator*) + vector arithmetic

namespace qchem::BasisSet::Lattice_3D
{

rvec3_t PeriodicGridEvaluator::GetGCartesian(const ivec3_t& m) const
{
    return itsRecip.GetCell().ToCartesian(rvec3_t(m)); // B m
}

// Uniform N1xN2xN3 grid of FRACTIONAL coordinates r=(i1/N1,i2/N2,i3/N3) -- the XC real-space grid
// (uniform weight Omega/prod(N); dG.r = 2 pi dm.r_frac makes the rho(G)<->rho(r) transform a plain DFT).
std::vector<rvec3_t> PeriodicGridEvaluator::UniformGrid(const ivec3_t& n) const
{
    std::vector<rvec3_t> g;
    g.reserve(static_cast<size_t>(n.x)*n.y*n.z);
    for (int i1=0;i1<n.x;i1++)
        for (int i2=0;i2<n.y;i2++)
            for (int i3=0;i3<n.z;i3++)
                g.push_back(rvec3_t(i1/double(n.x), i2/double(n.y), i3/double(n.z)));
    return g;
}

// Cartesian points of the FFT grid, raster order matching RhoOnGrid/ForwardFFT: r = A (i/N), A = direct cell.
// Cached: the N-point mesh + the 3x3 cell inversion are built once, not per DoFit.
const rvec3vec_t& PeriodicGridEvaluator::GridPoints() const
{
    if (itsGridPoints.size()==0)
    {
        std::vector<rvec3_t> frac=UniformGrid(itsN);
        UnitCell A=itsRecip.GetCell().MakeReciprocalCell();   // reciprocal of B == the direct cell
        itsGridPoints.resize(frac.size());
        for (size_t q=0;q<frac.size();q++) itsGridPoints[q]=A.ToCartesian(frac[q]);
    }
    return itsGridPoints;
}

// rho(r) on the FFT grid = inverse FFT of rho-tilde: rho(r_j) = Sum_dm rho-tilde(dm) e^{+i2pi dm.j/N}.
// rho-tilde is the physical coefficient (already /Omega), so the inverse FFT takes NO 1/N normalization.
rvec_t PeriodicGridEvaluator::RhoOnGrid(const ΔG_Map& rho) const
{
    ivec3_t N=itsN;
    size_t Npts=size_t(N.x)*N.y*N.z;
    cvec_t g(Npts, dcmplx(0.0));
    for (const auto& kv : rho)
    {
        const ivec3_t& m=kv.first;
        // TRUNCATE, don't ALIAS: a G beyond THIS grid's Nyquist (|m| > N/2) is not representable here, and the
        // modulo wrap below would FOLD it onto a low frequency (aliasing).  For an under-resolved sharp pair
        // (F's tight α≈80 density product, or a coarser XC grid receiving the collocator's multigrid top-rung
        // {G} that exceed its own Nyquist) that fold is a spurious CHECKERBOARD -> locally-negative ρ -> the XC
        // ρ>0 guard's grid-sensitive Exc collapse (doc/GPWPlan §0e step 2: NaF neg-frac 0.5, negCharge −9 e).
        // Dropping it is the graceful band-limit (CP2K-like); a grid-resolved density (Si: 0.08% negative) is
        // unaffected, and the per-level integrate-back only ever passes its OWN {G} (all in-band -> no drop).
        const int ax=m.x<0?-m.x:m.x, ay=m.y<0?-m.y:m.y, az=m.z<0?-m.z:m.z;
        if (2*ax>N.x || 2*ay>N.y || 2*az>N.z) continue;
        int i0=((m.x%N.x)+N.x)%N.x, i1=((m.y%N.y)+N.y)%N.y, i2=((m.z%N.z)+N.z)%N.z;
        g[(size_t(i0)*N.y+i1)*N.z+i2]=kv.second;
    }
    cvec_t rr=qchem::FFT::FFT3D(g, N, +1);
    rvec_t out(Npts);
    for (size_t i=0;i<Npts;i++) out[i]=std::real(dcmplx(rr[i]));
    return out;
}

// Forward-FFT a real-space grid field to the FULL normalised G-space grid Vtilde = FFT[V]/Npts (raster order).
// Unlike a difference-filtered map, this keeps EVERY resolved frequency, so any consumer's difference set hits.
cvec_t PeriodicGridEvaluator::ForwardFFT(const rvec_t& V) const
{
    ivec3_t N=itsN;
    size_t Npts=size_t(N.x)*N.y*N.z;
    assert(V.size()==Npts);
    cvec_t g(Npts, dcmplx(0.0));
    for (size_t i=0;i<Npts;i++) g[i]=dcmplx(V[i]);
    cvec_t Vt=qchem::FFT::FFT3D(g, N, -1);
    for (size_t i=0;i<Npts;i++) Vt[i]/=double(Npts);
    return Vt;
}

// Complex-input forward FFT (general-k): identical convention, the field is already complex so no real->complex
// pack.  FFT3D takes its argument by value, so pass a copy.
cvec_t PeriodicGridEvaluator::ForwardFFT(const cvec_t& V) const
{
    ivec3_t N=itsN;
    size_t Npts=size_t(N.x)*N.y*N.z;
    assert(V.size()==Npts);
    cvec_t Vt=qchem::FFT::FFT3D(V, N, -1);
    for (size_t i=0;i<Npts;i++) Vt[i]/=double(Npts);
    return Vt;
}

// Inverse of ForwardFFT: dense physical coefficients (raster order) -> real field, no normalisation
// (RhoOnGrid's FFT core without the sparse-map placement).
rvec_t PeriodicGridEvaluator::BackwardFFT(const cvec_t& c) const
{
    ivec3_t N=itsN;
    size_t Npts=size_t(N.x)*N.y*N.z;
    assert(c.size()==Npts);
    cvec_t rr=qchem::FFT::FFT3D(c, N, +1);
    rvec_t out(Npts);
    for (size_t i=0;i<Npts;i++) out[i]=std::real(dcmplx(rr[i]));
    return out;
}

// Vtilde(dm) from a ForwardFFT grid: wrap dm into [0,N) per axis (negative freq -> N-|.|) and index raster.
dcmplx PeriodicGridEvaluator::GridCoeff(const cvec_t& Vt, const ivec3_t& dm) const
{
    ivec3_t N=itsN;
    // Alias-free ONLY while |dm| < N/2 per axis: outside that band the wrap folds a high frequency onto dm
    // and silently returns the wrong coefficient.  Holds by construction for relCutoff>=1 (denser fit grid
    // is protective); the real future exposure is a k!=0 XC block whose orbital difference set exceeds the
    // Gamma fit grid -- this assert makes that path fail loudly instead of aliasing.
    assert(2*(dm.x<0?-dm.x:dm.x)<N.x && 2*(dm.y<0?-dm.y:dm.y)<N.y && 2*(dm.z<0?-dm.z:dm.z)<N.z
           && "GridCoeff: |dm| outside the alias-free band -- densify the fit grid (raise relCutoff)");
    int i0=((dm.x%N.x)+N.x)%N.x, i1=((dm.y%N.y)+N.y)%N.y, i2=((dm.z%N.z)+N.z)%N.z;
    return Vt[(size_t(i0)*N.y+i1)*N.z+i2];
}

// integral f d3r on the FFT grid: uniform quadrature, weight Omega/Npts.
double PeriodicGridEvaluator::Integral(const rvec_t& f) const
{
    return blazem::sum(f)*itsVolume/double(f.size());
}

// Inverse-transform the Hermitian coefficients c(dm) to the real field f(r) = Re Sum_dm c(dm) e^{i(B.dm).r}
// using this engine's B (GetGCartesian).  A point evaluation (loops the sparse map), NOT an FFT.
double PeriodicGridEvaluator::EvalField(const ΔG_Map& c, const rvec3_t& r) const
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
rvec3_t PeriodicGridEvaluator::EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const
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

} //namespace
