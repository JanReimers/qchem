// doc/Records/enzyme_loopcache_repro.cpp -- Enzyme reverse mode returns a WRONG gradient, with no diagnostic, when a callee fills a
// loop-local stack array through an inner loop whose trip count is read from data.  Self-contained: no Enzyme header, no qchem code.
// Lane B, 2026-11-08; the record and the draft issue are in doc/Records/EnzymeSpike.md §6.
//
//   CXX=/opt/LLVM-21.1.6-Linux-X64/bin/clang++; P=$HOME/Code/Enzyme/build/Enzyme/ClangEnzyme-21.so
//   $CXX -std=c++17 -O2 -fplugin=$P enzyme_loopcache_repro.cpp && ./a.out
//   expected:  dE/dD_00: AD 4.5709930...  FD 4.5709930...
//   actual:    dE/dD_00: AD 0.8510516320  FD 4.5709930507        (the same at -O1, -O2, -O3, with or without -march=native)
// Replacing the inner k loop by its one iteration (or a compile-time trip count) gives the right answer.
#include <cmath>
#include <cstdio>
extern int enzyme_dup, enzyme_const;
template <typename R, typename... T> R __enzyme_autodiff(T...);

struct Data { int n, np; const int* nk; const double* x; };   // nk[i]: the inner loop's trip count for entry i

static void Fill(const Data* d, int p, double* phi)
{
    for (int i=0;i<d->n;i++)
    {
        double v=0.0;
        for (int k=0;k<d->nk[i];k++) v+=std::exp(-(0.5+i+k)*d->x[p]);   // trip count read from data
        phi[i]=v;
    }
}
double E(const double* D, const Data* d)   // E(D) = sum_p (phi_p^T D phi_p)^2, phi_p independent of D
{
    double s=0.0;
    for (int p=0;p<d->np;p++)
    {
        double phi[3];
        Fill(d, p, phi);
        double rho=0.0;
        for (int i=0;i<d->n;i++) for (int j=0;j<d->n;j++) rho+=phi[i]*D[i*d->n+j]*phi[j];
        s+=rho*rho;
    }
    return s;
}
int main()
{
    const int nk[3]={1,1,1}; const double x[4]={0.1,0.7,1.3,2.0};
    const Data d{3, 4, nk, x};
    double D[9]={0.9,0.2,0.1, 0.2,0.7,0.0, 0.1,0.0,0.4}, g[9]={0};
    __enzyme_autodiff<void>((void*)E, enzyme_dup, D, g, enzyme_const, &d);
    const double h=1e-6; double Dp[9], Dm[9]; for (int k=0;k<9;k++) Dp[k]=Dm[k]=D[k]; Dp[0]+=h; Dm[0]-=h;
    std::printf("dE/dD_00: AD %.10f  FD %.10f\n", g[0], (E(Dp,&d)-E(Dm,&d))/(2*h));
}
