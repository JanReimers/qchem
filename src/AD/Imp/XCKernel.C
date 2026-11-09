// File: AD/Imp/XCKernel.C  The delta XC kernel and its Enzyme gradients (the formula: src/AD/XCKernel.C).
//
// Written for Enzyme, which differentiates the IR of THIS TU (the probes, 2026-11-08):
//  - every function the energy calls is defined here: a call into another TU or shared library is refused at compile time ("No
//    augmented forward pass found"), and so is a virtual call.  So the Becke cutoff and the functional are restated here, each beside
//    the library home it must equal (a unit test holds them equal), not imported.
//  - OpenMP stays OUTSIDE the differentiated code.  E is a sum over independent blocks of (atom, point) pairs, so the energy runs the
//    blocks in parallel and sums them in block order (deterministic), and a gradient differentiates ONE SERIAL BLOCK per call, the
//    calls running in parallel, their gradients summed in block order.  Enzyme's own OpenMP support failed on this kernel, three ways
//    (2026-11-08): it refuses a `reduction` clause (__kmpc_reduce_nowait), a vector allocated in the region segfaulted, and with the
//    partials passed in it corrupted the heap ("malloc(): unaligned tcache chunk") whenever the region captured a constant pointer
//    beside an active one (its "Mismatched activity" warning).  A simple loop with no constant captures did work (the spike's probe).
//  - a mesh node of zero template weight (MHL's r = 0 node) is skipped BEFORE any distance is taken: |r - R_owner| = 0 there, and the
//    derivative of a square root at zero, times a zero weight, is a NaN in the reverse pass.
module;
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include <enzyme/enzyme>
module qchem.AD.XCKernel;

namespace qchem::AD
{
namespace
{
constexpr size_t kMaxAtoms=64, kMaxBF=512;   // the per-point stack tables (checked once, at the entry points)

//! qchem.Mesh's BeckeCutoff (Mesh.C), restated: s(mu) = (1 - f^{k+1}(mu))/2, f(mu) = (3 mu - mu^3)/2.
inline double Cutoff(double mu, int k)
{
    for (int i=0;i<=k;i++) mu=0.5*(3.0*mu-mu*mu*mu);
    return 0.5*(1.0-mu);
}
//! qchem.Mesh.BeckePartition's BeckePartitionWeight, restated on plain arrays: atom \a owner's weight at the point \a r.
inline double Becke(const double* R, size_t n, size_t owner, const double* r, int k)
{
    if (n==1) return 1.0;
    double dist[kMaxAtoms];
    for (size_t i=0;i<n;i++)
    {
        const double dx=r[0]-R[3*i], dy=r[1]-R[3*i+1], dz=r[2]-R[3*i+2];
        dist[i]=std::sqrt(dx*dx+dy*dy+dz*dz);
    }
    double sum=0.0, own=0.0;
    for (size_t i=0;i<n;i++)
    {
        double P=1.0;
        for (size_t j=0;j<n;j++) if (j!=i)
        {
            const double dx=R[3*i]-R[3*j], dy=R[3*i+1]-R[3*j+1], dz=R[3*i+2]-R[3*j+2];
            const double Rij=std::sqrt(dx*dx+dy*dy+dz*dz);
            const double mu = Rij>0.0 ? (dist[i]-dist[j])/Rij : 0.0;
            P*=Cutoff(mu,k);
        }
        sum+=P;
        if (i==owner) own=P;
    }
    return sum>0.0 ? own/sum : 0.0;
}

//! The unpolarized LDA energy density per particle, eps_x + eps_c: SlaterExchange(2/3) (eps_x = 3/4 v_x, v_x = -3 alpha (3 rho_s/4 pi)^(1/3)
//! with rho_s = rho/2, in its operation order) + VWN_Correlation's paramagnetic G(x; P), x = sqrt(r_s).
constexpr double kPi=3.14159265358979323846, kFourPi=4.0*kPi;
inline double EpsXc(double rho)
{
    const double vx=-3.0*(2.0/3.0)*std::pow(3.0*(0.5*rho)/kFourPi, 1.0/3.0);
    const double A=0.0310907, b=3.72744, c=12.9352, x0=-0.10498;
    const double x=std::sqrt(std::cbrt(3.0/(kFourPi*rho)));
    const double Q=std::sqrt(4.0*c-b*b), Xx=x*x+b*x+c, X0=x0*x0+b*x0+c;
    const double at=std::atan(Q/(2.0*x+b)), beta=b*x0/X0;
    const double t1=std::log(x*x/Xx)             + (2.0*b/Q)*at;
    const double t2=std::log((x-x0)*(x-x0)/Xx)   + (2.0*(b+2.0*x0)/Q)*at;
    return 0.75*vx + A*(t1-beta*t2);
}

inline double IntPow(double x, int n) { double p=1.0; for (int i=0;i<n;i++) p*=x; return p; }

//! The basis values at \a r: Cartesian Gaussians riding with their atoms.
inline void AO(const XCKernelData& d, const double* R, const double* r, double* phi)
{
    for (const XCShell& s : d.shells)
    {
        const double dx=r[0]-R[3*s.atom], dy=r[1]-R[3*s.atom+1], dz=r[2]-R[3*s.atom+2], r2=dx*dx+dy*dy+dz*dz;
        double rad=0.0;
        for (size_t k=0;k<s.exponents.size();k++) rad+=s.coefficients[k]*std::exp(-s.exponents[k]*r2);
        for (size_t c=0;c<s.powers.size();c++)
            phi[s.offset+c]=s.norms[c]*IntPow(dx,s.powers[c][0])*IntPow(dy,s.powers[c][1])*IntPow(dz,s.powers[c][2])*rad;
    }
}

//! ONE BLOCK of the energy, serial: the (atom, point) pairs [npair b / nb, npair (b+1) / nb).  What a gradient differentiates.
double BlockEnergy(const double* R, const double* D, const XCKernelData& d, int b)
{
    const size_t nq=d.w.size(), npair=d.natom*nq, nbf=d.nbf, nb=size_t(d.nBlocks);
    const size_t lo=npair*size_t(b)/nb, hi=npair*size_t(b+1)/nb;
    double s=0.0;
    for (size_t p=lo;p<hi;p++)
    {
        const size_t A=p/nq, q=p%nq;
        if (d.w[q]==0.0) continue;                  // MHL's r = 0 node: before any distance (the file header)
        const double r[3]={R[3*A]+d.u[3*q], R[3*A+1]+d.u[3*q+1], R[3*A+2]+d.u[3*q+2]};
        const double wB=Becke(R, d.natom, A, r, d.beckeOrder);
        if (wB==0.0) continue;
        double phi[kMaxBF];
        AO(d, R, r, phi);
        double rho=0.0;
        for (size_t i=0;i<nbf;i++)
        {
            double t=0.0;
            for (size_t j=0;j<nbf;j++) t+=D[i*nbf+j]*phi[j];
            rho+=phi[i]*t;
        }
        if (rho>0.0) s+=d.w[q]*wB*rho*EpsXc(rho);
    }
    return s;
}

//! THE KERNEL: E_xc(R, D), the blocks in parallel, summed in block order.
double Energy(const double* R, const double* D, const XCKernelData& d)
{
    std::vector<double> part(d.nBlocks, 0.0);
#pragma omp parallel for schedule(dynamic)
    for (int b=0;b<d.nBlocks;b++) part[b]=BlockEnergy(R, D, d, b);
    double E=0.0;
    for (double e : part) E+=e;
    return E;
}

} // namespace

// The kernels Enzyme differentiates, and their constant contexts.  NOT in the anonymous namespace: Enzyme's entry point is a function
// template DECLARED and never defined (the plugin supplies it), and a template used with an argument type of internal linkage is
// ill-formed ("used but not defined ... its type does not have linkage").  Module linkage is enough.
namespace xc_detail
{
// The two block kernels: one block's energy, active in R or in D.
struct Ctx { const XCKernelData* d; const double* R; const double* D; int block; };
double KernelR(const double* R, const Ctx* c) {return BlockEnergy(R, c->D, *c->d, c->block);}
double KernelD(const double* D, const Ctx* c) {return BlockEnergy(c->R, D, *c->d, c->block);}

struct BeckeCtx { size_t n, owner; double u[3]; int k; };
double KernelBecke(const double* R, const BeckeCtx* c)
{
    const double r[3]={R[3*c->owner]+c->u[0], R[3*c->owner+1]+c->u[1], R[3*c->owner+2]+c->u[2]};
    return Becke(R, c->n, c->owner, r, c->k);
}
} // namespace xc_detail

namespace
{
using namespace xc_detail;

//! THE TYPED WRAPPER: the reverse-mode gradient of \a F(x, c) in the \a n doubles at \a x, \a c held constant, ACCUMULATED into \a g.
//! \a F is a template argument because Enzyme needs the function as a compile-time constant.
template <auto F, class C>
void ReverseInto(const double* x, double* g, const C* c)
{
    // Enzyme's CLASSIC entry: every argument's activity stated explicitly (enzyme_dup / enzyme_const).  The untyped function pointer
    // is Enzyme's ABI and stops here.  (Its tuple wrapper enzyme::autodiff is equivalent: a zero gradient once blamed on it was Enzyme's
    // loop-cache defect, doc/Records/EnzymeSpike.md §4.)
    __enzyme_autodiff<void>(reinterpret_cast<void*>(F), enzyme_dup, x, g, enzyme_const, c);
}

//! The gradient of the whole energy in the \a n doubles at \a x: one block per call, the calls in parallel (each into its own row),
//! the rows summed in block order.  \a C gets the block index from \a ctx(b).
template <auto F, class MakeCtx>
std::vector<double> BlockGradient(const double* x, size_t n, int nBlocks, MakeCtx ctx)
{
    std::vector<double> rows(size_t(nBlocks)*n, 0.0);
#pragma omp parallel for schedule(dynamic)
    for (int b=0;b<nBlocks;b++) { const auto c=ctx(b); ReverseInto<F>(x, &rows[size_t(b)*n], &c); }
    std::vector<double> g(n, 0.0);
    for (int b=0;b<nBlocks;b++) for (size_t i=0;i<n;i++) g[i]+=rows[size_t(b)*n+i];
    return g;
}

void Check(const XCKernelData& d, const std::vector<double>& R)
{
    if (d.natom>kMaxAtoms) throw std::invalid_argument("XCKernel: "+std::to_string(d.natom)+" atoms; the kernel's tables hold "+std::to_string(kMaxAtoms));
    if (d.nbf>kMaxBF)      throw std::invalid_argument("XCKernel: "+std::to_string(d.nbf)+" basis functions; the kernel's tables hold "+std::to_string(kMaxBF));
    if (R.size()!=3*d.natom || d.D.size()!=d.nbf*d.nbf || d.u.size()!=3*d.w.size() || d.nBlocks<1)
        throw std::invalid_argument("XCKernel: inconsistent data (R, D, mesh template or block count)");
    for (const XCShell& s : d.shells)
        if (s.atom>=d.natom || s.offset+s.powers.size()>d.nbf || s.norms.size()!=s.powers.size() || s.coefficients.size()!=s.exponents.size())
            throw std::invalid_argument("XCKernel: a shell does not fit the basis");
}
} // namespace

double XCEnergy(const XCKernelData& d, const std::vector<double>& R)
{
    Check(d,R);
    return Energy(R.data(), d.D.data(), d);
}

rvec_t XCEnergyGradientR(const XCKernelData& d, const std::vector<double>& R)
{
    Check(d,R);
    // No alias of the active R in the constant context (Enzyme's activity analysis).
    const auto g=BlockGradient<KernelR>(R.data(), R.size(), d.nBlocks, [&](int b){ return Ctx{&d, nullptr, d.D.data(), b}; });
    rvec_t out(R.size());
    for (size_t i=0;i<R.size();i++) out[i]=g[i];
    return out;
}

rmat_t XCEnergyGradientD(const XCKernelData& d, const std::vector<double>& R)
{
    Check(d,R);
    const std::vector<double> D=d.D;   // a COPY: the active D is not the constant context's d.D
    const auto g=BlockGradient<KernelD>(D.data(), D.size(), d.nBlocks, [&](int b){ return Ctx{&d, R.data(), nullptr, b}; });
    rmat_t out(d.nbf, d.nbf);
    for (size_t i=0;i<d.nbf;i++) for (size_t j=0;j<d.nbf;j++) out(i,j)=g[i*d.nbf+j];
    return out;
}

std::vector<double> XCKernelAOValues(const XCKernelData& d, const std::vector<double>& R, const rvec3_t& r)
{
    Check(d,R);
    std::vector<double> phi(d.nbf, 0.0);
    const double p[3]={r.x, r.y, r.z};
    AO(d, R.data(), p, phi.data());
    return phi;
}

double XCKernelBeckeWeight(const std::vector<double>& R, size_t owner, const rvec3_t& u, int k)
{
    const size_t n=R.size()/3;
    if (n>kMaxAtoms || owner>=n) throw std::invalid_argument("XCKernelBeckeWeight: owner or atom count out of range");
    const BeckeCtx c{n, owner, {u.x, u.y, u.z}, k};
    return KernelBecke(R.data(), &c);
}

rvec_t XCKernelBeckeWeightGradient(const std::vector<double>& R, size_t owner, const rvec3_t& u, int k)
{
    const size_t n=R.size()/3;
    if (n>kMaxAtoms || owner>=n) throw std::invalid_argument("XCKernelBeckeWeightGradient: owner or atom count out of range");
    const BeckeCtx c{n, owner, {u.x, u.y, u.z}, k};
    std::vector<double> g(R.size(), 0.0);
    ReverseInto<KernelBecke>(R.data(), g.data(), &c);
    rvec_t out(R.size());
    for (size_t i=0;i<R.size();i++) out[i]=g[i];
    return out;
}

} // namespace qchem::AD
