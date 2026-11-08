// doc/Records/enzyme_flags_repro.cpp -- Enzyme 1a87ccc (built against LLVM 21.1.6) differentiating the AD spike's XC block kernel
// (src/AD/Imp/XCKernel.C as of 2026-11-08, restated without modules), dE/dD by reverse mode against a central finite difference.
// Lane B, 2026-11-08; the record is doc/Records/EnzymeSpike.md.
//
//   P=$HOME/Code/Enzyme/build/Enzyme/ClangEnzyme-21.so; I=$HOME/Code/Enzyme/enzyme/include; CXX=/opt/LLVM-21.1.6-Linux-X64/bin/clang++
//   $CXX -std=c++20 -O2                                      -fplugin=$P -I$I enzyme_flags_repro.cpp && ./a.out   # AD == FD
//   $CXX -std=c++20 -O3                                      -fplugin=$P -I$I enzyme_flags_repro.cpp && ./a.out   # AD WRONG, silently
//   $CXX -std=c++20 -O2 -march=native                        -fplugin=$P -I$I enzyme_flags_repro.cpp && ./a.out   # SEGFAULT (AVX2 box)
//   $CXX -std=c++20 -O2 -march=native -mprefer-vector-width=128 -fplugin=$P -I$I enzyme_flags_repro.cpp && ./a.out   # AD == FD
// The line printed is "full <k> AD <value> FD <value>" for three entries of D.
#include <enzyme/enzyme>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
struct XCShell { size_t atom=0, offset=0; std::vector<std::array<int,3>> powers; std::vector<double> norms, exponents, coefficients; };
struct XCKernelData { size_t natom=0, nbf=0; int beckeOrder=3; std::vector<double> u, w; std::vector<XCShell> shells; std::vector<double> D; int nBlocks=256; };
namespace K {
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

//! THE KERNEL: E_xc(R, D).  \a part: d.nBlocks doubles of scratch, one partial sum per block (the file header).
double Energy(const double* R, const double* D, const XCKernelData& d, double* part)
{
    const size_t nq=d.w.size(), npair=d.natom*nq, nbf=d.nbf;
    const int nb=d.nBlocks;
//SERIALTEST #pragma omp parallel for schedule(static)
    for (int b=0;b<nb;b++)
    {
        const size_t lo=npair*size_t(b)/size_t(nb), hi=npair*size_t(b+1)/size_t(nb);
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
        part[b]=s;
    }
    double E=0.0;
    for (int b=0;b<nb;b++) E+=part[b];
    return E;
}

}

// The kernels Enzyme differentiates, and their constant contexts.  NOT in the anonymous namespace: Enzyme's entry point is a function
// template DECLARED and never defined (the plugin supplies it), and a template used with an argument type of internal linkage is
// ill-formed ("used but not defined ... its type does not have linkage").  Module linkage is enough.
namespace xc_detail
{
using namespace K;
// The two energy kernels: the same energy, active in R or in D.
struct Ctx { const XCKernelData* d; const double* R; const double* D; };
double KernelR(const double* R, const Ctx* c, double* part) {return Energy(R, c->D, *c->d, part);}
double KernelD(const double* D, const Ctx* c, double* part) {return Energy(c->R, D, *c->d, part);}

struct BeckeCtx { size_t n, owner; double u[3]; int k; };
double KernelBecke(const double* R, const BeckeCtx* c, double*)
{
    const double r[3]={R[3*c->owner]+c->u[0], R[3*c->owner+1]+c->u[1], R[3*c->owner+2]+c->u[2]};
    return Becke(R, c->n, c->owner, r, c->k);
}
} // namespace xc_detail
using namespace K; using namespace xc_detail;
int main(){
    XCKernelData d; d.natom=2; d.beckeOrder=3; d.nBlocks=3;
    const double u[][3]={{0.3,0.1,-0.2},{-0.7,0.9,0.4},{1.1,-0.5,0.8},{0.2,0.2,0.2},{-0.4,-0.6,0.1}};
    for (auto& p : u) { d.u.push_back(p[0]); d.u.push_back(p[1]); d.u.push_back(p[2]); d.w.push_back(0.37); }
    d.shells.push_back({0, 0, {{0,0,0}}, {1.3}, {0.8, 2.5}, {0.6, 0.4}});
    d.shells.push_back({1, 1, {{0,0,0}}, {0.9}, {1.1}, {1.0}});
    d.shells.push_back({0, 2, {{1,0,0},{0,1,0},{0,0,1}}, {1.7,1.7,1.7}, {0.9}, {1.0}});
    d.nbf=5;
    const double Dm[5][5]={{0.9,0.2,0.1,0.0,0.05},{0.2,0.7,0.0,0.1,0.0},{0.1,0.0,0.4,0.0,0.0},{0.0,0.1,0.0,0.3,0.02},{0.05,0.0,0.0,0.02,0.2}};
    for (auto& r : Dm) for (double v : r) d.D.push_back(v);
    std::vector<double> R={0.0,0.0,0.0, 0.4,0.3,1.6};
    Ctx c{&d, R.data(), nullptr};
    std::vector<double> D=d.D, g(25,0.0), part(3), dpart(3,0.0);
    __enzyme_autodiff<void>((void*)KernelD, enzyme_dup, D.data(), g.data(), enzyme_const, &c, enzyme_dup, part.data(), dpart.data());
    double h=1e-6; for(int k : {0,1,6}){ auto p=D,m=D; p[k]+=h; m[k]-=h; std::vector<double> pp(3),pm(3); std::printf("full %d AD %.10f FD %.10f\n",k,g[k],(KernelD(p.data(),&c,pp.data())-KernelD(m.data(),&c,pm.data()))/(2*h)); }
}
