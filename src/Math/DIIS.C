// File: DIIS.C  Pulay's DIIS bordered-system solve — the shared extrapolator engine.
//
// DIIS (Direct Inversion of the Iterative Subspace): P. Pulay, Chem. Phys. Lett. 73, 393 (1980);
// P. Pulay, J. Comput. Chem. 3, 556 (1982).  This is the ONE paper-faithful engine that serves BOTH
// Fock-space DIIS (extrapolate F from a history of [F,D] errors) and density-space Pulay/Broyden mixing
// (extrapolate ρ̃ from a history of ρ_out−ρ_in residuals).  Only the *residual stream* differs (its vectors
// and inner product); the linear algebra here is identical — see doc/SCFStrategyPlan.md §4.
//
// Given the error-overlap matrix B (Bᵢⱼ = ⟨eᵢ,eⱼ⟩, symmetric, n×n) of a history of n residual vectors,
// find the coefficients c (Σ cᵢ = 1) that minimise ‖Σ cᵢ eᵢ‖².  The Lagrange condition is the BORDERED
// system (Pulay 1980, Eq (6)):
//     ⎡ B   1 ⎤ ⎡ c ⎤   ⎡ 0 ⎤
//     ⎣ 1ᵀ  0 ⎦ ⎣ λ ⎦ = ⎣ 1 ⎦      (B_bordered is (n+1)×(n+1): B top-left, a border of 1, a 0 corner)
// Solve, then drop λ.  The caller owns the history + builds B with ITS OWN inner product.
module;
export module qchem.Math.DIIS;
import qchem.Blaze;   // rsmat_t / rmat_t / rvec_t + blazem::{zero,svd,solve,subvector}
import qchem.Math;

export namespace qchem::Math::DIIS
{

//! Build the (n+1) BORDERED matrix from the raw n×n error-overlap B (Pulay 1980, Eq (6)).
inline rsmat_t Bordered(const rsmat_t& B)
{
    size_t n=B.rows();
    rsmat_t Bb=blazem::zero<double>(n+1);
    for (size_t i=0;i<n;++i)
    {
        Bb(i,n)=1.0;                             // the border of 1 (symmetric storage mirrors Bb(n,i))
        for (size_t j=i;j<n;++j) Bb(i,j)=B(i,j); // B in the top-left (upper triangle; rsmat_t is symmetric)
    }
    return Bb;                                   // Bb(n,n)=0 already
}

//! Smallest singular value of the bordered system — the DIIS conditioning gauge (prune history when it drops).
inline double MinSV(const rsmat_t& Bb)
{
    rvec_t s; rmat_t U,Vt;
    blazem::svd(Bb,U,s,Vt);
    return s[s.size()-1];
}

//! Solve the bordered system for the coefficients c (Pulay 1980, Eq (6)).  RHS = e_last (the Σcᵢ=1 row).
inline rvec_t Coefficients(const rsmat_t& Bb)
{
    size_t N=Bb.rows();                          // n+1
    rvec_t rhs(N,0.0); rhs[N-1]=1.0;             // [0,…,0,1]
    rvec_t c=blazem::solve(Bb,rhs);
    return blazem::subvector(c,0,N-1);           // the n coefficients c (drop λ)
}

} //namespace
