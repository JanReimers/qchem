# Battery Materials — Critical-Path Roadmap

The connective tissue above the individual plans ([[doc/SCFSeedingPlan]],
[[doc/MolecularPseudopotentialPlan]], [[doc/GaussianPlaneWavePlan]]). North-star:
Li/Na cathode voltage curves (LiCoO₂ / LiMn₂O₄ / LiNiO₂) — [[project_battery_voltage_goal]].

**Detail decreases outward on purpose.** Near stages are specced; far stages are sketched, because
implementing the near ones will reshape the far ones. The value of the far sketches is *forward
visibility* — §"Forward-design hooks" lets the near-term implementer design better interfaces today
for stages two steps out, without over-specifying them.

## The reframe: two tiers (USPP is the upgrade, not the gate)

You can get a **first voltage curve with the GTH norm-conserving PPs you already have**, provided the
*physics* is adequate for transition-metal oxides. USPP/PAW is the efficiency+accuracy upgrade, not a
prerequisite. So correctness physics (GGA, +U) comes *before* the fancy pseudopotential.

- **Tier 1 — first voltage curve (minimum viable):** solid DFT (GPW) + GGA + spin + DFT+U + forces +
  cluster-expansion/MC, on NC PPs. Expensive (high cutoff) but *correct*.
- **Tier 2 — VASP-class:** USPP/PAW (soft pseudization → affordable) + k-point BZ throughput.

## The pipeline (dependency-ordered; detail decreasing outward)

```
0. IN FLIGHT / PLANNED   qcMath · SCF-seeding(SAD) · Molecular-PPs · GPW
1. PHYSICS ADEQUACY      GGA(∇ρ) · spin-in-solids · DFT+U              ← first voltages possible here
2. EFFICIENT PP          augmentation facet → USPP → (PAW)            ← the "VASP/USPP" target
3. FORCES                Hellmann-Feynman + Pulay (incl. dV_PP/dR)    ← interleaved with 1-2
4. THROUGHPUT            k-point sampling + Tier-A space-group BZ reduction
5. STAT-MECH             cluster expansion → Monte Carlo → V(T,x)     ← the actual deliverable
```

### Stage 0 — basis/Coulomb infrastructure (planned)
SAD seeding (ionic convergence), molecular PPs (the agnostic-interface forcing function), GPW
(Gaussian orbitals + PW Coulomb; the molecule↔solid unifier). See the three plans.

### Stage 1 — physics adequacy (the "right vs runs" stage)
- **GGA (∇ρ):** libxc has the functionals; the work is producing `∇ρ` on the grid and the GGA potential
  (carries a `∇·` term).
- **Spin in solids:** the Pol axis exists; make it work in the periodic/GPW path (magnetic Co/Ni/Mn).
- **DFT+U:** essential for localized 3d (plain LDA/GGA gets voltages badly wrong). A Hamiltonian term
  acting on **atom-centered d-orbital occupations** — reuses the localized-projector machinery (see hooks).

### Stage 2 — efficient pseudopotentials (the named target)
- **Augmentation** is the new concept: `S = 1 + Σ Q|β⟩⟨β|` (overlap stops being trivial) and the density
  gains a compensation charge `ρ = ρ_smooth + ρ_aug`. **USPP** = NC-separable **+ an augmentation facet**;
  **PAW** adds the all-electron reconstruction.
- **Structural break worth knowing:** `S ≠ I` is **not new machinery for the Gaussian/GPW path** —
  Gaussians already solve a generalized eigenproblem every SCF. Only the *pure-PW* path assumed `S = I`.
  So **GPW + USPP is more natural than PW + USPP** — GPW is on the battery path, not a detour.

### Stage 3 — forces
Hellmann-Feynman + Pulay, including `dV_PP/dR`. Needed to relax each Li configuration before it enters
the cluster-expansion training set. Realistically interleaved with Stages 1–2.

### Stage 4 — throughput: k-points + Tier-A space groups
- **k-point BZ sampling** for accurate energetics / small-gap systems.
- **Tier-A space-group reduction** (from the qcSymmetry TODO): use the crystal **point group acting on k
  + time-reversal (k ↔ −k)** to fold the k-mesh to the **irreducible wedge** with weights → far fewer
  k-points per run. This is the throughput multiplier for the many supercell DFT runs the CE training set
  needs. Builds on the molecular point-group/SALC core ([[project_molecular_symmetry_plan]]), extended to
  k-space. (**Tier B** — full `{R|τ}` non-symmorphic space groups, crystal SALCs, density symmetrization —
  is *not* on the battery path; defer as an accuracy/elegance item.)

### Stage 5 — statistical mechanics (the voltage curve)
Cluster expansion: fit a lattice-gas (Ising 0/1) Hamiltonian to the DFT total energies of many Li
configs. Monte Carlo at `(T, μ)` → `V(T, x)`. A separate module downstream of converged TEs + forces.

## Forward-design hooks (the point of writing the far stages now)

**Molecular-PP stage (Stage 0) should design for:**
- **USPP (2):** do NOT bake `S = I` or `ρ = Σ D χχ` into the PP capability or density assembly. USPP is
  NC-separable **+ an augmentation facet** (Q charges → `S` and `ρ`); the ISP/facet structure already
  takes "another facet" cleanly — just don't foreclose it.
- **Shared projectors:** the separable **β-projectors are reused three ways** — NC-separable, USPP, and
  **DFT+U's d-orbital projection** (Stage 1). Keep the projector interface a general atom-centered
  localized projector, not GTH-bound.
- **Forces (3):** keep `dV_PP/dR` computable — mesh-quadrature local → grid derivative; separable →
  `dβ/dR`. Don't strand the derivative path.

**GPW stage (Stage 0) should design for:**
- **GGA (1):** collocation must be able to emit `{ρ, ∇ρ}`, not just `ρ`.
- **USPP (2):** the grid density carries an additive augmentation charge `ρ_grid = ρ_smooth + ρ_aug`.
- **Forces (3):** the collocation adjoint supports position derivatives.

**k-sampling (whenever general-k lands, GPW or PW) should design for:**
- **Tier-A (4):** carry **per-k weights from the start** and make the BZ sum weight-respecting, so
  symmetry reduction folds in without restructuring the k-loop.

## Critical-path summary
First voltage curve is gated by **GGA + spin + DFT+U + forces + CE/MC on NC PPs** — *not* by USPP.
USPP/PAW + Tier-A k-reduction are the Tier-2 upgrades that take it from "works" to "VASP-class
throughput." Molecular-PPs and GPW (Stage 0) are the infrastructure both tiers stand on, which is why
their interfaces should already see USPP, +U, GGA, and forces coming.
