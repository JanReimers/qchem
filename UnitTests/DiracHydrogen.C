// Dirac hydrogen atom with RKB Gaussian basis (even-tempered)
#include <blaze/Math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

import qchem.LASolver;
import qchem.Types;
import qchem.Blaze;

using std::cout;
using std::endl;

static double gaussian_norm(double alpha) {
    return std::pow(2.0*alpha/M_PI, 0.75);
}

static double overlap_ij(double ai, double aj) {
    double p = ai + aj;
    return std::pow(M_PI / p, 1.5);
}
static double nuc_attraction_ij(double ai, double aj, double Z) {
    double p = ai + aj;
    return -2.0 * M_PI * Z / p;
}
static double deriv_ij(double ai, double aj) {
    double p = ai + aj;
    return -4.0 * M_PI * aj / (p*p);
}
static double invr_ij(double ai, double aj) {
    double p = ai + aj;
    return 2.0 * M_PI / p;
}
static double kinetic_ij(double ai, double aj) {
    // primitive kinetic integral for s-type GTOs: <g_i| -1/2 nabla^2 |g_j>
    double p = ai + aj;
    return 3.0 * ai * aj * std::pow(M_PI, 1.5) / (2.0 * std::pow(p, 2.5));
}

void dump(const rsmat_t& S, size_t N, const char* name, size_t ioff=0,size_t joff=0)
{
    std::cout << name << "=" << endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << S(ioff+i, joff+j) << " ";
        std::cout << "\n";
    }    
}

int main() {
    const double Z = 1.0;
    const double c_light = 137.035999084;   // atomic units
    const double m = 1.0;
    const double mc2 = m * c_light * c_light;

    const double alpha0 = 0.01;
    const double beta = 2.0;
    const int N = 22;

    LASolver<double>* las=LASolver<double>::Factory(qchem::Cholsky,1e-12);

    std::vector<double> alphas(N);
    for (int i = 0; i < N; ++i) alphas[i] = alpha0 * std::pow(beta, i);

    std::vector<double> norm(N);
    for (int i = 0; i < N; ++i) norm[i] = gaussian_norm(alphas[i]);

    rsmat_t S=zero<double>(N), V=zero<double>(N),T_nonrel=zero<double>(N);

    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            double Nij = norm[i]*norm[j];
            S(i,j) = Nij * overlap_ij(alphas[i], alphas[j]);
            V(i,j) = Nij * nuc_attraction_ij(alphas[i], alphas[j], Z);
            // non-relativistic kinetic integral <phi| -1/2 nabla^2 |phi'>
            T_nonrel(i,j) = Nij * kinetic_ij(alphas[i], alphas[j]);
        }
    }

    const int M = 2 * N;
    rsmat_t H=zero<double>(M), Sbig=zero<double>(M);
    rsmat_t Tmat=zero<double>(M), Vmat=zero<double>(M), Rmat=zero<double>(M);

    // Build RKB small-component via p^2 = 2 * T_nonrel
    rsmat_t P = 2.0 * T_nonrel; // P_{ij} = <phi_i|p^2|phi_j>
    const double pref_ss = 1.0 / (4.0 * m * m * c_light * c_light);
    const double pref_ls = 1.0 / (2.0 * m);

    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            double Nij = norm[i]*norm[j];
            // Large-large
            H(i, j) = V(i, j);
            Sbig(i, j) = S(i, j);
            Vmat(i, j) = V(i, j);
            Rmat(i, j) =0.0;
            // Small-small using RKB: S_SS = (1/4m^2 c^2) * P
            Sbig(N+i, N+j) = pref_ss * P(i,j);
            // compute V_SS = (1/(4 m^2 c^2)) * <phi_i| p V p |phi_j>
            // Using RKBS small-component 1/r matrix element formula:
            //   Inv_r1 = 4*alpha_i*alpha_j*Integral(t,2*l+1)*ns_i*ns_j  (already includes 4*pi)
            // For l=0, Integral(t,1) = 1/(2*p), so Vss = -Z * 2*alpha_i*alpha_j/p * ns_i*ns_j * pref_ss.
            double p = alphas[i] + alphas[j];
            double Vss_prim = -2.0 * Z * alphas[i] * alphas[j] / p;
            double Vss = pref_ss * Nij * Vss_prim;
            Vmat(N+i, N+j) = Vss;
            
            Rmat(N+i, N+j) = -2*mc2 * Sbig(N+i, N+j);
            H(i, j) = V(i, j) + Rmat(i, j);
            H(N+i, N+j) = Vmat(N+i, N+j) + Rmat(N+i, N+j);
            // RKB coupling between large and small component bases
            Tmat(i, N+j) = pref_ls * P(i,j);
            Tmat(N+i, j) = pref_ls * P(i,j);
            H(i, N+j) = Tmat(i, N+j);
            H(N+i, j) = Tmat(N+i, j);
        }
    }

    //dump small-component blocks for debug comparison
    const int dumpN = std::min(3, N);
    std::cout << "\nDebug dump: first " << dumpN << " elements\n";
    dump(Sbig,dumpN,"S_LL",0,0);
    dump(Sbig,dumpN,"S_SS",N,N);
    dump(H   ,dumpN,"H_LS",0,N);
    dump(H   ,dumpN,"H_SL",N,0);
    dump(Vmat,dumpN,"V_LL",0,0);
    dump(Vmat,dumpN,"V_SS",N,N);

    las->SetBasisOverlap(Sbig);
    auto [U,w]=las->Solve(H); //U is rmat_t, w is rvec_t. U is already back transformed.
    
    if ((int)w.size() != M) {
        std::cerr << "Eigen decomposition failed\n";
        return 3;
    }
    // cout << "All eigen values=" << w << endl;

    // collect real eigenpairs and map back to original coefficients
    struct EV { double val; int idx; };
    std::vector<EV> electrons;
    for (int i = 0; i < M; ++i) {
        // discard positron solutions with E < -mc2
        if (w[i] > -mc2) electrons.push_back({w[i], i});
    }
    // cout << "Electron eigen values={";
    // for (auto e:electrons) cout << e.val << ",";
    // cout << "}" << endl;

    if (electrons.empty()) {
        std::cerr << "No electron eigenvalues found after filtering positrons\n";
        return 4;
    }

    // pick ground-state electron as the one with smallest (E - mc2)
    double bestVal = 1e300; int bestIdx = -1;
    for (auto &e: electrons) {
        double Eb = e.val - mc2;
        if (Eb < bestVal) { bestVal = Eb; bestIdx = e.idx; }
    }

    rvec_t y = column(U, bestIdx);
    double E_y = std::real(w[bestIdx]);
    
    // bestIdx corresponds to eigenvalue index in w
    double E_eig_full = w[bestIdx];
    double E_bind  = E_eig_full - mc2;

    double Ekin  = trans(y) * Tmat * y;
    double Epot  = trans(y) * Vmat * y;
    double Erest = trans(y) * Rmat * y;
    double E_H   = trans(y) * H * y;
    double Etot = Ekin + Epot + Erest;
    
    double max_resid = blaze::max(blaze::abs(H * y - E_eig_full * y));
    
    double Etot_bind = Etot - mc2;

    std::cout.setf(std::ios::fixed); std::cout.precision(8);
    std::cout << "Found electron states: " << electrons.size() << " / total " << M << "\n";
    std::cout << "Selected eigenvalue E = " << E_eig_full << "  (E-mc2 = " << E_bind << ")\n";
    std::cout << "Expectation values: T = " << Ekin << ", V = " << Epot << ", R = " << Erest << ", T+V+R = " << Etot << "\n";
    std::cout << "Direct H expectation = " << E_H << "\n";
    std::cout << "Residual max |Hc - E S c| = " << max_resid << "\n";
    std::cout << "Expectation (Ebinding) = " << Etot_bind << "  reference = -0.50000665\n";
    std::cout << "Difference eigen-binding - reference = " << (E_bind - (-0.50000665)) << "\n";
    std::cout << "Difference expectation-binding - reference = " << (Etot_bind - (-0.50000665)) << "\n";

    return 0;
}
