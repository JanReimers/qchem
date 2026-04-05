// dirac_hydrogen_gauss.cpp
#include <blaze/Math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>

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

int main() {
    const double Z = 1.0;
    const double c_light = 137.035999084;   // au
    const double m = 1.0;
    const double mc2 = m * c_light * c_light;
    const int N = 5;
    const std::vector<double> alphas = {0.25, 0.75, 2.5, 8.0, 25.0};

    if ((int)alphas.size() != N) {
        std::cerr << "N mismatch\n"; return 2;
    }

    std::vector<double> norm(N);
    for (int i = 0; i < N; ++i) norm[i] = gaussian_norm(alphas[i]);

    blaze::DynamicMatrix<double> S(N,N,0.0), V(N,N,0.0);
    blaze::DynamicMatrix<double> Kplus(N,N,0.0), Kminus(N,N,0.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double Nij = norm[i]*norm[j];
            double Sij = overlap_ij(alphas[i], alphas[j]);
            double Vij = nuc_attraction_ij(alphas[i], alphas[j], Z);
            S(i,j) = Nij * Sij;
            V(i,j) = Nij * Vij;

            double d = deriv_ij(alphas[i], alphas[j]);
            double invr = invr_ij(alphas[i], alphas[j]);
            Kplus(i,j)  = Nij * (d - invr);
            Kminus(i,j) = Nij * (-d - invr);
        }
    }

    const int M = 2 * N;
    blaze::DynamicMatrix<double> H(M,M,0.0), Sbig(M,M,0.0);
    blaze::DynamicMatrix<double> Tmat(M,M,0.0), Vmat(M,M,0.0);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            H(i, j) = V(i, j) + mc2 * S(i, j);
            H(N+i, N+j) = V(i, j) - mc2 * S(i, j);
            Sbig(i, j) = S(i, j);
            Sbig(N+i, N+j) = S(i, j);
            Vmat(i, j) = V(i, j);
            Vmat(N+i, N+j) = V(i, j);

            Tmat(i,j) = mc2 * S(i,j);
            Tmat(N+i,N+j) = -mc2 * S(i,j);
            Tmat(i, N+j) = c_light * Kminus(i, j);
            Tmat(N+i, j) = c_light * Kplus(i, j);
        }
    }

    // Enforce exact symmetry
    H = 0.5 * (H + trans(H));
    Sbig = 0.5 * (Sbig + trans(Sbig));
    Tmat = 0.5 * (Tmat + trans(Tmat));
    Vmat = 0.5 * (Vmat + trans(Vmat));

    // Check symmetry
    double dH = max(abs(H - trans(H)));
    double dS = max(abs(Sbig - trans(Sbig)));
    std::cout << "H symmetry error = " << dH << "\n";
    std::cout << "S symmetry error = " << dS << "\n";

    blaze::DynamicMatrix<double> S_inv = blaze::inv(Sbig);
    blaze::DynamicMatrix<double> A = S_inv * H;
    blaze::DynamicMatrix<double> A_sym = 0.5 * (A + trans(A));

    blaze::DynamicVector<std::complex<double>> w;
    blaze::DynamicMatrix<std::complex<double>> U;
    blaze::eigen(A_sym, w, U);

    if (w.size() != (size_t)M) {
        std::cerr << "Eigen failure\n"; return 3;
    }

    std::vector<double> realEv;
    for (size_t i = 0; i < w.size(); ++i) {
        double im = std::abs(w[i].imag());
        if (im > 1e-9) {
            std::cerr << "Warning: small imaginary part in eigenvalue: " << w[i] << "\n";
        }
        realEv.push_back(w[i].real());
    }
    std::sort(realEv.begin(), realEv.end());

    std::cout << "Eigenvalues:\n";
    for (int i = 0; i < M; ++i) {
        double E = realEv[i];
        std::cout << "  " << i << " E = " << E << " (E-mc2 = " << (E-mc2) << ")\n";
    }

    // filter out positron states.
    // ground-state eigenvector index (lowest energy)
    int gidx = 0;
    for (int i = 1; i < M; ++i)
        if (realEv[i] < realEv[gidx] && realEv[i] > 0) gidx = i;

    // find in generalized eigen decomposition
    int idxU = -1;
    for (int i = 0; i < M; ++i) {
        if (std::abs(w[i].real() - realEv[gidx]) < 1e-8) { idxU = i; break; }
    }
    if (idxU < 0) { std::cerr << "No ground vector found\n"; return 4; }

    blaze::DynamicVector<std::complex<double>> c = column(U, idxU);
    blaze::DynamicVector<std::complex<double>> cn = c / std::sqrt( std::real( trans(c) * Sbig * c ) );
    double realNorm = std::real( trans(cn) * Sbig * cn );
    std::cout << "Normalized S-overlap of selected GS vector = " << realNorm << "\n";

    double Ekin = std::real( trans(cn) * Tmat * cn );
    double Epot = std::real( trans(cn) * Vmat * cn );
    double Etot = Ekin + Epot;

    std::cout << "Ground state diagonalization energy = " << realEv[gidx] << "\n";
    std::cout << "Expectation: T = " << Ekin << ", V = " << Epot << ", T+V = " << Etot << "\n";
    std::cout << "Rest-mass shifted: E-mc2 = " << (realEv[gidx] - mc2)
              << ", expectation T+V-mc2 = " << (Etot - mc2) << "\n";
    std::cout << "Benchmark E-mc2 (~1s Dirac H) = -0.49999 a.u.\n";

    return 0;
}