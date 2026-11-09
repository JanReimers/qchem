// Enzyme reverse mode SEGFAULTS with -march=native (AVX2).
//
// E(D) = sum over grid points of w * rho * eps(rho), rho = phi^T D phi, with phi a few Cartesian Gaussians on two atoms, a Becke
// partition weight per point, and an LDA energy density eps.  The point loop is split into 3 blocks; each block writes its partial
// sum into a scratch buffer `part` (passed enzyme_dup), and the partials are summed after the loop.  dE/dD by __enzyme_autodiff is
// compared with a central finite difference for three entries of D.
//
//   clang++ -std=c++20 -O2                                         -fplugin=ClangEnzyme-21.so -I<enzyme>/include repro.cpp && ./a.out
//       AD == FD
//   clang++ -std=c++20 -O2 -march=native                           -fplugin=ClangEnzyme-21.so -I<enzyme>/include repro.cpp && ./a.out
//       Segmentation fault (AVX2 machine)
//   clang++ -std=c++20 -O2 -march=native -mprefer-vector-width=128 -fplugin=ClangEnzyme-21.so -I<enzyme>/include repro.cpp && ./a.out
//       AD == FD
//   clang++ -std=c++20 -O3                                         -fplugin=ClangEnzyme-21.so -I<enzyme>/include repro.cpp && ./a.out
//       AD wrong, no diagnostic (AD -0.0130734452 against FD -1.1419468942 for D_00)
//
// With nBlocks = 1 in main(), -O2 -march=native is correct.  Smaller cases cut down from this file did not crash.
#include <enzyme/enzyme>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
struct Shell { size_t atom=0, offset=0; std::vector<std::array<int,3>> powers; std::vector<double> norms, exponents, coefficients; };
struct Data { size_t natom=0, nbf=0; int beckeOrder=3; std::vector<double> u, w; std::vector<Shell> shells; std::vector<double> D; int nBlocks=256; };
namespace K {
constexpr size_t kMaxAtoms=64, kMaxBF=512;

// Becke's cutoff: s(mu) = (1 - f^{k+1}(mu))/2, f(mu) = (3 mu - mu^3)/2.
inline double Cutoff(double mu, int k)
{
    for (int i=0;i<=k;i++) mu=0.5*(3.0*mu-mu*mu*mu);
    return 0.5*(1.0-mu);
}
// Becke partition weight of atom `owner` at the point r.
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

// LDA energy density per particle: Slater exchange + VWN paramagnetic correlation.
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

// Basis function values at r: Cartesian Gaussians centred on the atoms.
inline void AO(const Data& d, const double* R, const double* r, double* phi)
{
    for (const Shell& s : d.shells)
    {
        const double dx=r[0]-R[3*s.atom], dy=r[1]-R[3*s.atom+1], dz=r[2]-R[3*s.atom+2], r2=dx*dx+dy*dy+dz*dz;
        double rad=0.0;
        for (size_t k=0;k<s.exponents.size();k++) rad+=s.coefficients[k]*std::exp(-s.exponents[k]*r2);
        for (size_t c=0;c<s.powers.size();c++)
            phi[s.offset+c]=s.norms[c]*IntPow(dx,s.powers[c][0])*IntPow(dy,s.powers[c][1])*IntPow(dz,s.powers[c][2])*rad;
    }
}

// E(R, D).  part: d.nBlocks doubles of scratch, one partial sum per block.
double Energy(const double* R, const double* D, const Data& d, double* part)
{
    const size_t nq=d.w.size(), npair=d.natom*nq, nbf=d.nbf;
    const int nb=d.nBlocks;
    for (int b=0;b<nb;b++)
    {
        const size_t lo=npair*size_t(b)/size_t(nb), hi=npair*size_t(b+1)/size_t(nb);
        double s=0.0;
        for (size_t p=lo;p<hi;p++)
        {
            const size_t A=p/nq, q=p%nq;
            if (d.w[q]==0.0) continue;
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

// Outside an anonymous namespace: __enzyme_autodiff is a template the plugin supplies, and its argument types need linkage.
namespace detail
{
using namespace K;
struct Ctx { const Data* d; const double* R; const double* D; };
double KernelD(const double* D, const Ctx* c, double* part) {return Energy(c->R, D, *c->d, part);}
} // namespace detail
using namespace K; using namespace detail;
int main(){
    Data d; d.natom=2; d.beckeOrder=3; d.nBlocks=3;
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
    double h=1e-6; for(int k : {0,1,6}){ auto p=D,m=D; p[k]+=h; m[k]-=h; std::vector<double> pp(3),pm(3); std::printf("D_%d: AD %.10f FD %.10f\n",k,g[k],(KernelD(p.data(),&c,pp.data())-KernelD(m.data(),&c,pm.data()))/(2*h)); }
}
