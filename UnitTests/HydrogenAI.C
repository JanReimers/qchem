// dirac_hydrogen_gauss.cpp
#include <blaze/Math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

static double gaussian_norm(double alpha) {
    // s-type GTO in 3D: N = (2*alpha/pi)^(3/4)
    return std::pow(2.0*alpha/M_PI, 0.75);
}

static double overlap_ij(double alpha_i, double alpha_j) {
    double p = alpha_i + alpha_j;
    return std::pow(M_PI / p, 1.5);
}

static double kinetic_ij(double alpha_i, double alpha_j) {
    double p = alpha_i + alpha_j;
    return 3.0 * alpha_i * alpha_j * std::pow(M_PI, 1.5) / (2.0 * std::pow(p, 2.5));
}

static double nuc_attraction_ij(double alpha_i, double alpha_j, double Z) {
    double p = alpha_i + alpha_j;
    return -2.0 * M_PI * Z / p;
}

//  ∫ φ_i (d/dr) φ_j d^3r for s-type radial envelope (3D)
static double deriv_ij(double alpha_i, double alpha_j) {
    double p = alpha_i + alpha_j;
    return -4.0 * M_PI * alpha_j / (p*p);
}

//  ∫ φ_i (1/r) φ_j d^3r for s-type GTO in same center
static double invr_ij(double alpha_i, double alpha_j) {
    double p = alpha_i + alpha_j;
    return 2.0 * M_PI / p;
}

int main() {
    const double Z = 1.0;
    const double c = 137.035999084;  // atomic units
    const double m = 1.0;
    const double mc2 = m * c * c;
    const int N = 5;  // number of GTO primitives for large component

    // Exponents chosen for 1s-like coverage
    const std::vector<double> alphas = {0.25, 0.75, 2.5, 8.0, 25.0};

    if ((int)alphas.size() != N) {
        std::cerr << "Change N to match alphas size.\n";
        return 1;
    }

    // Precompute normalization
    std::vector<double> norm(N);
    for (int i = 0; i < N; ++i)
        norm[i] = gaussian_norm(alphas[i]);

    // Build non-relativistic blocks
    blaze::DynamicMatrix<double> S(N, N, 0.0);  // overlap
    blaze::DynamicMatrix<double> V(N, N, 0.0);  // nuclear
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double Sij = overlap_ij(alphas[i], alphas[j]);
            double Nij = norm[i] * norm[j];
            S(i,j) = Nij * Sij;
            V(i,j) = Nij * nuc_attraction_ij(alphas[i], alphas[j], Z);
        }
    }

    // Dirac off-diagonal coupling (L-S and S-L)
    blaze::DynamicMatrix<double> Kplus(N, N, 0.0);   // d/dr - 1/r
    blaze::DynamicMatrix<double> Kminus(N, N, 0.0);  // -(d/dr + 1/r)

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double Nij = norm[i] * norm[j];
            double d = deriv_ij(alphas[i], alphas[j]);
            double invr = invr_ij(alphas[i], alphas[j]);
            // kappa = -1 for 1s (j=1/2, l=0)
            // O_plus = d/dr - 1/r
            Kplus(i,j) = Nij * (d - invr);
            // O_minus = -(d/dr + 1/r) = -d/dr - 1/r
            Kminus(i,j) = Nij * (-d - invr);
        }
    }

    // Build Dirac 2N x 2N Hamiltonian and overlap
    const int M = 2 * N;
    blaze::DynamicMatrix<double> H(M, M, 0.0);
    blaze::DynamicMatrix<double> Sbig(M, M, 0.0);

    // LL and SS diagonal blocks
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            H(i, j) = V(i, j) + mc2 * S(i, j);             // H_LL
            H(N+i, N+j) = V(i, j) - mc2 * S(i, j);         // H_SS
            Sbig(i, j) = S(i, j);
            Sbig(N+i, N+j) = S(i, j);
        }
    }

    // Off-diagonal couplings
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            H(i, N+j) = c * Kminus(i, j);
            H(N+i, j) = c * Kplus(i, j);
        }
    }

    {
        blaze::DynamicMatrix<double> U,Vt;
        blaze::DynamicVector<double> S_svd;
        blaze::svd(Sbig, U, S_svd, Vt);
        std::cout << "Smallest singular value of S: " << S_svd[2*N-1] << "\n";
    }
    // Generalized eigen problem H c = E Sbig c
    // Convert to standard eigen by computing Sinv * H
    blaze::DynamicMatrix<double> S_inv = blaze::inv(Sbig);
    blaze::DynamicMatrix<double> A = S_inv * H;

    std::cout << "Assymetry in H=" << max(abs(H-trans(H))) << std::endl;
    std::cout << "Assymetry in A=" << max(abs(A-trans(A))) << std::endl;

    // Solve (generally non-Hermitian due approx), but eigenvalues should be real
    blaze::DynamicVector<std::complex<double>> w;
    blaze::DynamicMatrix<std::complex<double>> U;
    blaze::eigen(A, w, U);

    // Pick lowest positive-energy root near mc2, and lowest below mc2 for binding
    std::vector<double> reals;
    for (size_t i = 0; i < w.size(); ++i)
        if (std::abs(w[i].imag()) < 1e-8)
            reals.push_back(w[i].real());
        else
            std::cout << "Complex eigenvalue: " << w[i] << "\n";

    std::sort(reals.begin(), reals.end());

    // Print spectrum
    std::cout << "Dirac spectrum (2N=" << M << " states):\n";
    for (size_t idx = 0; idx < reals.size(); ++idx) {
        double E = reals[idx];
        double eps = E - mc2;  // binding relative to rest energy
        std::cout << "  #" << (idx+1)
                  << "  E = " << E
                  << "  (E-mc2 = " << eps << " a.u.)\n";
    }
    std::cout.flush();
    if (reals.empty()) {
        std::cerr << "No real eigenvalues found!\n";
        return 1;
    }
    // In good Dirac 1s: E - mc2 ≈ -0.5 a.u. (Z=1)
    double best = reals.front();
    if (!reals.empty()) {
        double binding = best - mc2;
        std::cout << "\nEstimated 1s energy (lowest): E = " << best
                  << ", E-mc2 = " << binding << " a.u.\n";
        std::cout << "Benchmark Dirac 1s (Z=1) approx E-mc2 = -0.5 a.u. (nonrel) and with relativity ~ -0.499992\n";
    }

    return 0;
}