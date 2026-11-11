// File: Structure/Ewald.C  Ewald summation of the point-ion electrostatic energy of a periodic cell.
//
// Computes the electrostatic energy of a set of point charges {qₐ at τₐ} repeated on the lattice of a
// UnitCell, immersed in a uniform neutralising background of charge −Q (Q=Σqₐ).  This is the ion–ion
// (nuclear–nuclear) term that a plane-wave total energy needs to become a physical (negative) cohesive
// energy: the bare lattice sum Σ' qₐqᵦ/|τₐᵦ+R| is only conditionally convergent, so it is split (Ewald)
// into a real-space sum (erfc-screened, short range) plus a reciprocal-space sum (the smooth part, a
// Fourier series over G), a per-ion self-energy correction, and a background term for the net charge.
//
// The split parameter η is a pure convergence knob: the total is η-INDEPENDENT (a regression/physics
// invariant the tests pin).  All quantities are atomic units (e=1, 4πε₀=1); energy in Hartree.
//
// The charges are the ION CORE charges — always positive (Si: 4,4; NaF: 1,7) — NOT formal oxidation
// states.  Covalent vs ionic bonding lives entirely in the SCF electron density, never in this term;
// so the SAME routine serves any crystal (covalent Si or ionic NaF), differing only in the data passed.
module;
#include <cassert>
#include <functional>
#include <vector>
#include <algorithm>
export module qchem.Ewald;
import qchem.Types;     // rvec_t, rvec3_t, vec3_t
import qchem.UnitCell;
import qchem.LatticePoints;   // LatticePointsInSphere (V-STRUCTURE-ISP F2: images off the cell)

namespace qchem {

//! Electrostatic (Madelung) energy per cell of the point charges \a q sitting at the UnitCell's atom
//! positions, lattice-repeated in a uniform neutralising background.  \a q is indexed like the cell's
//! atoms (q[a] is the charge of cell[a]); pass the ION CORE charges (valence/pseudo-ion charge).
//! \a eta is the Ewald split (a.u.⁻¹); eta<=0 picks a balanced default.  The result is eta-independent.
export double EwaldEnergy(const UnitCell& cell, const rvec_t& q, double eta=0.0)
{
    const size_t N = cell.GetNumAtoms();
    assert(q.size()==N);
    const double Ω = cell.GetCellVolume();

    // Cartesian basis positions and the charge moments.
    std::vector<rvec3_t> τ(N);
    double Q=0.0, Qsq=0.0;
    for (size_t a=0; a<N; a++)
    {
        τ[a]  = cell.Sites()[a].R;
        Q    += q[a];
        Qsq  += q[a]*q[a];
    }

    // Balanced default split: η ~ 1/Ω^{1/3} keeps both the real and reciprocal sums small.
    if (eta<=0.0) eta = sqrt(Pi)/cbrt(Ω);

    // Cutoffs for a target accuracy e^{-p}: erfc(η·rcut)~e^{-p}, exp(-gcut²/4η²)~e^{-p}.
    const double p    = 34.0;             // ~e^{-34} ≈ 2e-15 truncation (η-independent to ~1e-10)
    const double rcut = sqrt(p)/eta;
    const double gcut = 2.0*eta*sqrt(p);

    // --- Real-space sum: ½ Σ_R Σ_{a,b} qₐqᵦ erfc(η r)/r,  excluding the a=b,R=0 self term ----------
    // The PAIR distance r = |τₐ−τᵦ+R| is what the cutoff bounds, so R must run to rcut + max|τₐ−τᵦ| (doc/Pins.md pin 29): enumerating
    // |R| ≤ rcut alone dropped images with r < rcut whenever two home atoms sit far apart -- 47 mHa on MnO AFM-II, whose O at
    // (¾,¾,¾) is 21.8 bohr from the Mn at the origin against rcut 21.9 (OpenWork §3 "MnO accuracy", 2026-10-23).
    double span=0.0;
    for (size_t a=0; a<N; a++) for (size_t b=0; b<N; b++) span=std::max(span, norm(τ[a]-τ[b]));
    double Er=0.0;
    const UnitCell& dir = cell;
    for (const vec3_t<int>& n : LatticePointsInSphere(dir, rcut+span))
    {
        rvec3_t R = dir.ToCartesian(rvec3_t(n.x,n.y,n.z));
        bool R0 = (n.x==0 && n.y==0 && n.z==0);
        for (size_t a=0; a<N; a++)
            for (size_t b=0; b<N; b++)
            {
                if (R0 && a==b) continue;          // drop the self term
                double r = norm(τ[a]-τ[b]+R);
                if (r<1e-12 || r>rcut) continue;   // coincident images guard; beyond the cutoff erfc is below e^{-p}
                Er += q[a]*q[b]*erfc(eta*r)/r;
            }
    }
    Er *= 0.5;

    // --- Reciprocal-space sum: (2π/Ω) Σ_{G≠0} e^{−G²/4η²}/G² |S(G)|²,  S(G)=Σₐ qₐ e^{iG·τₐ} ---------
    double Eg=0.0;
    const CellGeometry recip=cell.MakeReciprocal();   // columns of B = reciprocal lattice vectors
    for (const vec3_t<int>& n : LatticePointsInSphere(recip, gcut))
    {
        if (n.x==0 && n.y==0 && n.z==0) continue;       // G=0 handled by the background term
        rvec3_t G  = recip.ToCartesian(rvec3_t(n.x,n.y,n.z));
        double  G2 = G*G;
        double Sre=0.0, Sim=0.0;
        for (size_t a=0; a<N; a++)
        {
            double φ = G*τ[a];
            Sre += q[a]*cos(φ);
            Sim += q[a]*sin(φ);
        }
        Eg += exp(-G2/(4.0*eta*eta))/G2 * (Sre*Sre+Sim*Sim);
    }
    Eg *= 2.0*Pi/Ω;

    // --- Self-energy of the erfc screening, and the net-charge / neutralising-background term --------
    double Eself = -eta/sqrt(Pi) * Qsq;
    double Ebg   = -Pi/(2.0*Ω*eta*eta) * Q*Q;

    return Er + Eg + Eself + Ebg;
}

//! THE FORCE of \c EwaldEnergy on each ion, \f$F_a=-\partial E/\partial\tau_a\f$ (GeometryOptimizationPlan §13.2), with the SAME split, cutoffs
//! and image enumeration, so it is the derivative of the number the energy returns.  The self and background terms do not depend on the
//! positions.  Per term:
//!   real:        F_a += Σ_{b,R}' q_a q_b [erfc(ηr)/r² + (2η/√π) e^{-η²r²}/r] d/r,   d = τ_a − τ_b + R
//!   reciprocal:  F_a += (4π/Ω) q_a Σ_{G≠0} (e^{-G²/4η²}/G²) G [Re S sin(G·τ_a) − Im S cos(G·τ_a)]
//! The result is η-independent and sums to zero (a uniform translation moves nothing).
export std::vector<rvec3_t> EwaldForces(const UnitCell& cell, const rvec_t& q, double eta=0.0)
{
    const size_t N = cell.GetNumAtoms();
    assert(q.size()==N);
    const double Ω = cell.GetCellVolume();
    std::vector<rvec3_t> τ(N);
    for (size_t a=0; a<N; a++) τ[a]=cell.Sites()[a].R;
    if (eta<=0.0) eta = sqrt(Pi)/cbrt(Ω);
    const double p    = 34.0;             // the energy's truncation, kept identical
    const double rcut = sqrt(p)/eta;
    const double gcut = 2.0*eta*sqrt(p);
    std::vector<rvec3_t> F(N, rvec3_t(0,0,0));

    double span=0.0;
    for (size_t a=0; a<N; a++) for (size_t b=0; b<N; b++) span=std::max(span, norm(τ[a]-τ[b]));
    for (const vec3_t<int>& n : LatticePointsInSphere(cell, rcut+span))
    {
        const rvec3_t R = cell.ToCartesian(rvec3_t(n.x,n.y,n.z));
        const bool R0 = (n.x==0 && n.y==0 && n.z==0);
        for (size_t a=0; a<N; a++)
            for (size_t b=0; b<N; b++)
            {
                if (R0 && a==b) continue;
                const rvec3_t d = τ[a]-τ[b]+R;
                const double  r = norm(d);
                if (r<1e-12 || r>rcut) continue;
                const double s = q[a]*q[b]*(erfc(eta*r)/(r*r) + 2.0*eta/sqrt(Pi)*exp(-eta*eta*r*r)/r)/r;
                F[a] += s*d;
            }
    }

    const CellGeometry recip=cell.MakeReciprocal();
    for (const vec3_t<int>& n : LatticePointsInSphere(recip, gcut))
    {
        if (n.x==0 && n.y==0 && n.z==0) continue;
        const rvec3_t G  = recip.ToCartesian(rvec3_t(n.x,n.y,n.z));
        const double  G2 = G*G;
        double Sre=0.0, Sim=0.0;
        for (size_t a=0; a<N; a++) { const double φ=G*τ[a]; Sre+=q[a]*cos(φ); Sim+=q[a]*sin(φ); }
        const double w = 4.0*Pi/Ω*exp(-G2/(4.0*eta*eta))/G2;
        for (size_t a=0; a<N; a++)
        {
            const double φ=G*τ[a];
            F[a] += (w*q[a]*(Sre*sin(φ)-Sim*cos(φ)))*G;
        }
    }
    return F;
}

//! THE VIRIAL of \c EwaldEnergy, \f$V_{ab}=\partial E/\partial\varepsilon_{ab}\f$ under the homogeneous strain \f$\tau\to(1+\varepsilon)\tau\f$,
//! \f$A\to(1+\varepsilon)A\f$ (doc/Records/CellHalf_Design.md §5; user ruling §6.2, 2026-11-10): a GRADIENT (the pullback's sign), in Ha,
//! so the stress is \f$\sigma=V/\Omega\f$.  The SAME split, cutoffs and image enumeration as the energy, at FIXED \f$\eta\f$ (the total is
//! \f$\eta\f$-independent, so the fixed-\f$\eta\f$ derivative is the derivative).  Per term, with \f$d=\tau_a-\tau_b+R\f$, \f$G\to(1+\varepsilon)^{-\top}G\f$:
//!   real:        ½ Σ' q_a q_b φ'(r) d_a d_b / r,       φ(r)=erfc(ηr)/r
//!   reciprocal:  −δ_ab E_g + (2π/Ω) Σ_{G≠0} |S(G)|² f'(G²) (−2 G_a G_b),   f(G²)=e^{−G²/4η²}/G²   (S is strain-invariant: G·τ is)
//!   self:        0;   background ∝ 1/Ω: −δ_ab E_bg.
//! Symmetric.  An index pair (a, b) is (row, column) of the returned matrix, 1-based (Matrix3D).
export rmat3d_t EwaldVirial(const UnitCell& cell, const rvec_t& q, double eta=0.0)
{
    const size_t N = cell.GetNumAtoms();
    assert(q.size()==N);
    const double Ω = cell.GetCellVolume();
    std::vector<rvec3_t> τ(N);
    double Q=0.0;
    for (size_t a=0; a<N; a++) { τ[a]=cell.Sites()[a].R; Q+=q[a]; }
    if (eta<=0.0) eta = sqrt(Pi)/cbrt(Ω);
    const double p    = 34.0;             // the energy's truncation, kept identical
    const double rcut = sqrt(p)/eta;
    const double gcut = 2.0*eta*sqrt(p);
    double V[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    auto C=[](const rvec3_t& v, int k){return k==0 ? v.x : k==1 ? v.y : v.z;};

    double span=0.0;
    for (size_t a=0; a<N; a++) for (size_t b=0; b<N; b++) span=std::max(span, norm(τ[a]-τ[b]));
    for (const vec3_t<int>& n : LatticePointsInSphere(cell, rcut+span))
    {
        const rvec3_t R = cell.ToCartesian(rvec3_t(n.x,n.y,n.z));
        const bool R0 = (n.x==0 && n.y==0 && n.z==0);
        for (size_t a=0; a<N; a++)
            for (size_t b=0; b<N; b++)
            {
                if (R0 && a==b) continue;
                const rvec3_t d = τ[a]-τ[b]+R;
                const double  r = norm(d);
                if (r<1e-12 || r>rcut) continue;
                const double dphi = -erfc(eta*r)/(r*r) - 2.0*eta/sqrt(Pi)*exp(-eta*eta*r*r)/r;   // φ'(r)
                const double s = 0.5*q[a]*q[b]*dphi/r;
                for (int i=0;i<3;i++) for (int j=0;j<3;j++) V[i][j]+=s*C(d,i)*C(d,j);
            }
    }

    const CellGeometry recip=cell.MakeReciprocal();
    double Eg=0.0;
    for (const vec3_t<int>& n : LatticePointsInSphere(recip, gcut))
    {
        if (n.x==0 && n.y==0 && n.z==0) continue;
        const rvec3_t G  = recip.ToCartesian(rvec3_t(n.x,n.y,n.z));
        const double  G2 = G*G;
        double Sre=0.0, Sim=0.0;
        for (size_t a=0; a<N; a++) { const double φ=G*τ[a]; Sre+=q[a]*cos(φ); Sim+=q[a]*sin(φ); }
        const double S2=Sre*Sre+Sim*Sim, f=exp(-G2/(4.0*eta*eta))/G2, df=f*(-1.0/(4.0*eta*eta)-1.0/G2);
        Eg += f*S2;
        const double s = 2.0*Pi/Ω*S2*df*(-2.0);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) V[i][j]+=s*C(G,i)*C(G,j);
    }
    Eg *= 2.0*Pi/Ω;
    const double Ebg = -Pi/(2.0*Ω*eta*eta) * Q*Q;
    for (int i=0;i<3;i++) V[i][i]-=Eg+Ebg;
    return rmat3d_t(V[0][0],V[0][1],V[0][2], V[1][0],V[1][1],V[1][2], V[2][0],V[2][1],V[2][2]);
}

} // namespace qchem