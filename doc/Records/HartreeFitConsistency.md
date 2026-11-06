# Hartree through a fit: when is H the gradient of E?

RECORD, 2026-11-06 (user and controller, from the TeF₆ fit-basis work; measurements by lane S).  Why the δ Hartree route is built the
way it is.  Open work stays in `doc/OpenWork.md` (the "(ij|δ)" row); the measurements are in `lane-s-ij-delta-measurement.md` and
`doc/Records/TeF6EpsilonReport.md`.

## The result

A fit of ρ gives a Hamiltonian that is the true gradient of its energy exactly when the implied four-index integral is symmetric under
ij ↔ kl.  That holds if the same object sits on both pairs:

    (ij|kl) = O_ij^T K O_kl        K symmetric          consistent: H = dE/dD
    (ij|kl) = C_ij^T G^-1 C_kl     Dunlap               consistent, and second order
    (ij|kl) = C_ij^T S^-1 O_kl     Coulomb one side,    NOT consistent
                                   overlap the other

Symbols below.  Forces need the consistent shapes.

## The pieces

For a δ fit basis on mesh points g with weights w_g, and a pair density P_ij(r) = χ_i(r)χ_j(r):

- **q**, the charges: q_g = w_g ρ(r_g).  These are the fit coefficients c_f.  Linear in D: q = O·D.
- **O**, the overlap-type vector: O_{g,ij} = w_g P_ij(r_g) = ⟨δ_g|ij⟩.
- **K**, the discrete Coulomb operator: v_g = Σ_h K_gh q_h.  `CoulombOperator::Potential` is K acting on q.  Coulomb's law is linear, so
  every method is a matrix, stored or not.  Direct sum: 1/|r_g − r_h| with a self term on the diagonal.  Radial Poisson: moments, radial
  solve and site sum, never stored.  Solids: FFT⁻¹ · 4π/G² · FFT.
- **A_ij(r)**, the analytic potential of a pair at a point: A_ij(r_g) = (ij|δ_g), a nuclear-attraction integral.  Finite: there is no
  singularity anywhere in (ij|δ).  The only singular object is the δ–δ element (δ_g|δ_h) on its diagonal, which a normal fit does not need
  and a Coulomb-metric fit does.
- **C, S, G** for a general fit basis f: C = (ij|f), S = ⟨f|f′⟩, G = (f|f′).

v_H itself is finite and smooth, two derivatives smoother than ρ (∇²v = −4πρ).  The singular things are the kernel and the nuclear
potential.  Fitting v is therefore reasonable, which is why numerical-grid and plane-wave codes do it.

## Two forms of one integral

The exact matrix is one double integral; doing the inner integral analytically gives two exact expressions:

    (a)  J_ij = ∫ ρ(r) A_ij(r) dr          ρ goes on the fit basis
    (b)  J_ij = ∫ P_ij(r) v(r) dr          v goes on the fit basis

On a mesh they differ:

    J^a_ij = Σ_g q_g A_ij(r_g)  =  Σ_f c_f (ij|δ_f)
    J^b_ij = Σ_g O_{g,ij} v_g

H_ij = H_ji holds for both.  What fails is one level up: ∂H_ij/∂D_kl ≠ ∂H_kl/∂D_ij.

Both give the same energy, E = ½Tr(DJ) = ½ Σ_g q_g v_g (measured, to every digit printed), but with exact v the gradient of that E is
neither form: dE/dD = J^c = ½(J^a + J^b).  J^a alone has no energy function at all, so its SCF would stop at a point not stationary in
E, and forces from E would be wrong.

## What makes one form consistent

**One-sided, analytic.**  J^a implies (ij|kl) = Σ_g A_ij(r_g) O_{g,kl}: a Coulomb integral on one pair, an overlap on the other.  The
same holds for any overlap-metric ("normal") fit, Gaussian included: c = S⁻¹O·D, H = C c, implied Cᵀ S⁻¹ O.  The first-order energy
and this mismatch are two faces of one choice.

**Dunlap.**  c = G⁻¹C·D gives Cᵀ G⁻¹ C.  Symmetric, and second order, because the fit minimizes the quantity the energy depends on.  For δ
functions G is the singular (δ_g|δ_h), so a Dunlap δ fit does not exist.

**Discrete operator.**  v = K q gives H = Oᵀ K q and E = ½ qᵀ K q, so ∂E/∂D = ½ Oᵀ(K + Kᵀ) q, which is H exactly when K = Kᵀ.  The
implied integral is O_ijᵀ K O_kl, the same O on both sides, and swapping ij ↔ kl transposes K.  First order, but consistent.  This is the
route on main (`GridDensityFitter` over a `CoulombOperator`).

**Solids never met this.**  GPW has the discrete-operator shape: collocate, 4π/G², gather with the adjoint.  A plane wave is an
eigenfunction of the Coulomb operator, so overlap and Coulomb fits coincide.

The trade: A_ij is exact, and that exactness breaks the pairing.  An approximate K between two identical O's stays consistent.

## Measured (TeF₆ PP, fixed D, max|ΔJ| against exact J; J^b here uses EXACT v)

| rung | J^a | J^b | J^c |
|---|---|---|---|
| 80/29 | 3.1e-5 | 4.6e-4 | 2.3e-4 |
| 150/47 | 2.7e-8 | 4.9e-7 | 2.5e-7 |

J^a is about 15× better than J^b, but J^c carries half of J^b's error, so the consistent analytic form buys about 2× at N_pts × N_pairs
integrals per build.  Why J^a beats J^b when v is the smoother field is unmeasured.  Guess: the sharp factor in J^b is a single core pair,
which max|ΔJ| finds, and in J^a it is the same ρ for every element.

## K's symmetry, measured (2026-11-06, lane S; gate `CoulombOperatorSymmetry`)

`RadialPoissonCoulomb` expands each source about its own site, so K_gh and K_hg can differ across sites.  |q′ᵀKq − qᵀKq′| / qᵀKq on TeF₆'s
geometry with model charges: about 4e-7 at 80/29, 5e-8 to 1e-7 at 150/47 (l ≥ 14), under 1e-8 at 300/47 (l ≥ 20).  The direct sum is symmetric
to 1e-9.  Scaled by E_H (347 Ha) that is about 1.4e-4, 3e-5 and 3e-6 Ha against ε·N = 7e-5, so inside ε at the geometry rung.  This is a
size, not a force bound: that stays open until forces exist.  The analytic element stays parked.

## Context: which codes do which (from memory, not checked)

Form (a) with a Gaussian auxiliary basis: Gaussian RI-J codes (Turbomole, ORCA, deMon), consistent through the Coulomb metric.  Form (b):
DMol³ and FHI-aims (multipole ρ, Poisson, H on the grid), SIESTA and CP2K's GPW (uniform grid), plane-wave codes.  The nearest relative of a
one-sided analytic δ form is pseudospectral / COSX exchange.
