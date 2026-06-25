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
#include <vector>
export module qchem.Ewald;
import qchem.Types;     // rvec_t, rvec3_t, vec3_t
import qchem.UnitCell;
import qchem.Structure;  // Atom (positions)
import qchem.Math;       // erfc, exp, sqrt, cos, sin, Pi, cbrt

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
        τ[a]  = cell[a]->itsR;
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
    double Er=0.0;
    const UnitCell& dir = cell;
    for (const vec3_t<int>& n : dir.CellsInSphere(rcut))
    {
        rvec3_t R = dir.ToCartesian(rvec3_t(n.x,n.y,n.z));
        bool R0 = (n.x==0 && n.y==0 && n.z==0);
        for (size_t a=0; a<N; a++)
            for (size_t b=0; b<N; b++)
            {
                if (R0 && a==b) continue;          // drop the self term
                double r = norm(τ[a]-τ[b]+R);
                if (r<1e-12) continue;             // coincident images guard
                Er += q[a]*q[b]*erfc(eta*r)/r;
            }
    }
    Er *= 0.5;

    // --- Reciprocal-space sum: (2π/Ω) Σ_{G≠0} e^{−G²/4η²}/G² |S(G)|²,  S(G)=Σₐ qₐ e^{iG·τₐ} ---------
    double Eg=0.0;
    const UnitCell recip = cell.MakeReciprocalCell();   // columns of B = reciprocal lattice vectors
    for (const vec3_t<int>& n : recip.CellsInSphere(gcut))
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
