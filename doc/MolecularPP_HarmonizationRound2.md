# Molecular ↔ Plane-Wave Pseudopotential Harmonization — Round 2 (the road to GPW)

**Status:** planning doc for the *next* session. Self-contained. Author-owned (like
`doc/MolecularPP_HarmonizationFindings.md`, "Round 1"); **not** the user's
`doc/MolecularPseudopotentialPlan.md`. Round 1 got molecular pseudopotentials working and harmonized the
*fitting/assembly* seams; Round 2 closes the last structural divergences and lays the track to **GPW**
(Gaussian orbitals on a periodic lattice), which is the north-star that makes solids/battery work possible.

---

## 0. Orientation (read this first — the doc is otherwise standalone)

Two pseudopotential (PP) code paths exist today and they must eventually share as much machinery as possible:

- **Molecular** (`src/Hamiltonian`, `double`, real-space, finite geometry): the external PP is assembled in the
  **term** — `PP_Local` (local part `⟨χᵢ|V_loc(r)|χⱼ⟩`) and `PP_NonLocal` (Kleinman–Bylander separable part
  `⟨χᵢ|β_p Yₗₘ⟩`), each quadratured on the geometry's own integration mesh.
- **Plane-wave / lattice** (`src/BasisSet/Lattice_3D` + `src/Hamiltonian/…/PWTerms.C`, `dcmplx`, reciprocal
  space, periodic): the PP is assembled in the **basis** — `PlaneWave_IBS` realizes
  `Integrals_Pseudo<dcmplx>` (structure factor `e^{-iG·τ}`, `1/Ω`, form factors), reached from the
  `PW_Pseudo` term by a `dynamic_cast` across abstract faces.

They already share a **spine** (this is what Round 1 established or confirmed):

| Shared object / algorithm | Role |
|---|---|
| `LocalPotential` / `SeparablePotential` (`src/Pseudopotential`) | the PP **models** (real + reciprocal views); both sides consume the same model |
| `Structure::CreateIntegrationMesh(mp)` (pure virtual) | each geometry owns its most efficient real-space mesh: `Atom`→radial×angular, `Molecule`→Becke, `UnitCell`→uniform grid. PW never calls it (owns its G-grid). |
| `Structure::SumFormFactors(f)` + `isFinite()` | the neutral G=0 alignment seam (finite → 0; periodic → folds in `1/Ω`), no `dynamic_cast<UnitCell>` |
| `NuclearRepulsion(st, zionOf)` | one ion-ion algorithm: pair-sum (finite) / Ewald (periodic); both `Vnn` and `PW_IonIon` delegate |
| `Fit_ABS` factory + `FunctionFitter_{Scalar,Density}<T>` | density/potential fitting; **both** sides now obtain their fitter through `bs->Create{CD,Vxc}FitBasisSet(...)` (Round 1 Item 2 — PW no longer assumes `orbital == fit`; `FourierFunctionFitter` retired) |

The map (colour-coded; **✓** = harmonized, **◐** = partly done, amber = the divergences this doc plans):

![Molecular vs plane-wave PP — shared spine + divergence points](diagrams/pp_molecular_vs_pw.svg)

**What Round 1 finished** (so we don't re-open it): `PseudoG0Energy` eliminated; `Integrals_Pseudo<T>` reduced
to two clean matrix methods (PW-only, `<dcmplx>`); both molecular PP terms made geometry-neutral on
`CreateIntegrationMesh`; multi-species routing; spin-native `Ham_PP` (`FittedVxcPol`/`FittedVcorrPol`, now pinned
by `A_PP.Si_PP_U.Polarized`); the whole PW DFT-fit path routed through the basis factory (Hartree + XC), with the
uniform `UnitCell` mesh + `L_PP` term-level validation in place. Full detail: `MolecularPP_HarmonizationFindings.md`
§6. **The "LASolver blocker" of Round 1 §5.1 was a red herring — near-singular basis conditioning, not a solver
bug; use a well-conditioned basis (Slater/High).**

---

## 1. The four remaining "incidental" divergences

"Incidental" = an artifact of how the two sides grew up, not a reflection of real physics (contrast the
**principled** divergence real-space `V(r)` vs reciprocal `ṽ(G²)`, which we *keep*). Each is a harmonization /
hoisting target. They are listed easy→hard; **(A) and (D) are ultimately the same problem, and GPW is what forces
them** (see §2).

### (B) G=0 alignment — trivially unifiable, do it opportunistically
- **Now:** the PP G=0 term is `0` for a finite structure and `(N/Ω)·Σ_a FormFactorG0(Z_a)` for a periodic one.
  `PW_Pseudo::GetEnergy` already branches on `!isFinite()` and reads the sum from `Structure::SumFormFactors`.
  The molecular side carries a hardcoded `0.0`.
- **Target:** when a single PP term serves all geometries (see (D)), it carries one expression
  `Ealign = isFinite() ? 0 : N·SumFormFactors(FormFactorG0)`. No new abstraction — the seam already exists.
- **Effort:** trivial (folds into (D)). **Risk:** none (bit-identical; finite value is provably 0).

### (C) `IonIon<T>` — collapse `Vnn` and `PW_IonIon` into one template
- **Now:** `Vnn` (`rsmat_t`, molecular) and `PW_IonIon` (`chmat_t`, periodic) are **the same energy-only term** —
  both add no matrix contribution and delegate to `NuclearRepulsion(st, zionOf)`. They differ only in scalar type.
- **Target:** one `IonIon<T>` term (`T ∈ {double, dcmplx}`), parameterized by the `zionOf` callback
  (identity for all-electron, `Z→Zion` map for PP), carrying the G=0 alignment of (B). Mirrors how the two PP
  terms already share their model + mesh.
- **Effort:** small (both terms already thin + delegating). **Risk:** low; per-side bit-identical anchors exist
  (`Enn = Zion_O·Zion_Si/R` etc.). Watch: the two live in different libs today — put `IonIon<T>` where both the
  `double` and `dcmplx` term bases are visible.

### (A) Assembly-in-term vs assembly-in-basis — the headline
- **Now:** molecular assembles in the **term** (`PP_Local`/`PP_NonLocal` + `CreateIntegrationMesh` + generic
  `qcMesh` quadrature); PW assembles in the **basis** (`Integrals_Pseudo<dcmplx>` via `dynamic_cast`).
- **Key insight (do not "fix" by writing a molecular `Integrals_Pseudo<double>`):** the pure PW path is
  *principled* — reciprocal-space assembly on its own G-grid is the efficient, correct thing, and it should stay.
  The molecular term path is *also* principled and is already geometry-neutral. They **converge only when a
  real-space (Gaussian) basis lives on a lattice** — that basis will quadrature its PP with the *molecular term
  path on a `UnitCell` mesh*, exactly like a finite `Molecule`. That is GPW. So (A) is **not** a refactor to do
  in isolation; it is **resolved by (D)/GPW**.
- **Effort:** none standalone — it dissolves under (D). **Risk:** the trap is "unify by giving molecular an
  `Integrals_Pseudo<double>`" — explicitly rejected: it would drag reciprocal-space assumptions into the molecular
  basis for zero benefit.

### (D) Full lattice-PP **SCF** — the real increment
- **Now:** the molecular PP terms already assemble correctly on a `UnitCell`'s uniform mesh — proven bit-identical
  by `L_PP` (`UnitTests/L_PP.C`: the same Si valence Gaussian basis + GTH PP gives the same `PP_Local`/
  `PP_NonLocal` matrices whether the atom is a finite `Molecule` or centred in a large `UnitCell`). But there is
  **no periodic SCF** with a real-space basis: two things are missing.
  1. **Periodic Gaussians.** Overlap / kinetic / Hartree of Gaussians on a lattice need Bloch sums
     (`Σ_R e^{ik·R} χ(r−R)`) and lattice-summed two-centre integrals (minimum-image or Ewald-style long-range).
     This is the substantive new numerics.
  2. **The facade must preserve the concrete geometry.** `qchem::Calculation`'s ctor does
     `itsStructure = std::make_shared<Molecule>(st)` — it **deep-copies any structure to a `Molecule`, stripping
     periodicity**. A lattice calculation must keep the `UnitCell` (and its `CreateIntegrationMesh` → uniform grid,
     and `Vnn`→Ewald). This is a small but load-bearing facade change and a **prerequisite for GPW**.
- **Target:** a real-space-basis SCF on a `UnitCell` reusing `PP_Local`/`PP_NonLocal`/`IonIon` unchanged.
- **Effort:** large (this is the GPW body). **Risk:** medium-high (new periodic-integral numerics); de-risk with
  the term-level `L_PP` bit-identity already in hand and the empty-lattice / cosine-V / bare-Coulomb PW anchors.

---

## 2. Strategic roadmap: symmetry → GPW

The user's stated next-steps, with the dependency structure made explicit. **GPW (§2.4) is the payoff and the
forcing function for divergences (A)+(D).** The symmetry work (§2.1–2.2) is valuable in its own right and
partially independent of GPW.

```
        (2.1) symmorphic space groups in qcSymmetry
                 │
         ┌───────┴────────┐
         ▼                ▼
 (2.2a) BZ reduction   (2.2b) SALC with plane waves
 (irreducible wedge)   (star-of-G symmetry blocking)
         │                │
         └───────┬────────┘   (both are *uses* of the space group; neither gates GPW)
                 ▼
        (2.3) PlaneWave_IBS "bad habits" review   ← do BEFORE GPW
                 ▼
        (2.4) GPW  ── resolves divergences (A) + (D); needs facade-preserves-UnitCell (§1.D.2)
```

### 2.1 Symmorphic space groups (foundation)
Extend `qcSymmetry` (which already has molecular point groups + `Lattice_3D`/Bloch machinery — note a qcSymmetry
folder/namespace reorg into `Atom`/`Molecule`/`Lattice_3D` recently landed) with **symmorphic** space groups =
point group ⋉ lattice translations (no screw axes / glide planes). New pieces: the space-group operations
(point op + lattice translation), the **star of k**, the **little co-group** of k, and its small (irreducible)
representations. *Symmorphic-first is the right scoping:* the little group's reps are just ordinary point-group
reps — **no fractional-translation phase / ray (projective) representations**, which is exactly the complexity
that non-symmorphic groups add. Defer non-symmorphic to a later round.

### 2.2a Brillouin-zone reduction (irreducible wedge) — *independent of SALC*
Use the symmorphic space group to fold the k-point mesh to the **irreducible Brillouin zone** and carry per-k
weights, so a periodic SCF only diagonalizes symmetry-distinct k-points. This is the **most broadly useful**
symmetry payoff for solids (it cuts SCF cost directly) and needs only §2.1 — not SALC, not GPW. Clean correctness
check: total energy from the reduced+weighted mesh equals the full-mesh energy.

### 2.2b SALC with plane waves — *independent of GPW*
Symmetry-adapt the `{G}` plane-wave basis: build symmetry-adapted combinations within each **star of G**, block-
diagonalizing the PW Hamiltonian. Mirrors the existing molecular SALC (`SymmetryAdapt`, `ShellRep`,
`OperationRep`) but on `{G}` instead of AO shells. Good, self-contained validation of §2.1 with a crisp check
(symmetry-blocked PW-SCF == unadapted PW-SCF). **Note:** GPW's *variational* basis is Gaussians (adapted by the
point group + Bloch), so SALC-PW is **not** a hard prerequisite for GPW — it's the cheaper way to exercise §2.1.

### 2.3 `PlaneWave_IBS` "bad habits" review — do before GPW
`PlaneWave_IBS` was coded early and accreted responsibilities that don't belong on a basis. Round 1 already
extracted **one** (fitting: "the basis should not do fits — `qcFitting` does, through a real independent fit basis
set," now honoured via `PlaneWaveFit_IBS` + the factory seam). GPW will lean hard on `PlaneWave_IBS` (it *is* the
reciprocal-space grid engine), so **audit it for the remaining cruft first** — a short, high-leverage cleanup
before building on top. Candidate smells to check (verify against current code — some may already be gone):
- SRP violations: does the basis still own things that are really term/Hamiltonian or fitting concerns
  (density→grid, ρ̃→Hartree/FFT-XC assembly currently kept concrete on it "for now")?
- `dynamic_cast`-reached capabilities that could be a clean abstract face (the CLAUDE.md cast policy).
- Anything assuming Γ-only or a single k where the lattice generalization needs the full k-set.
- Ownership/mesh/grid duplication now that the fitter owns its own grid.
Output: a ranked list + the cheap ones done, mirroring the Round-1 fitting extraction.

### 2.4 GPW (Gaussian And Plane Waves) — the payoff
**Method (CP2K / Lippert–Hutter):** orbitals in Gaussians (compact, good for core/valence); represent the
electron **density on a regular real-space grid**; FFT to G-space and solve Poisson there for Hartree
(`V_H(G) = 4π ρ(G)/G²`); evaluate XC on the grid; integrate the grid potential back against the Gaussians to form
the KS matrix. This is precisely where molecular and lattice PP **become one code path**:
- The KS matrix element `∫ φ_i V_grid φ_j` **is** `PP_Local`'s `qcMesh::WeightedOverlap(mesh, basis, V)` pattern
  generalized to an arbitrary grid potential `V` — resolving divergence (A).
- Running that on a `UnitCell` with periodic Gaussians is divergence (D).

**What we already have that GPW reuses:** the uniform `UnitCell` real-space mesh; `qchem.FFT`; the reciprocal-space
Poisson/Hartree machinery in `PW_Hartree` (and the Γ PW fit basis + factory seam from Round 1); geometry-neutral
PP terms; `Vnn`→Ewald on a `UnitCell`. **What's genuinely new:** (1) periodic Gaussian two-centre integrals
(Bloch/lattice sums), (2) the **collocate/integrate** pair (Gaussian density-matrix → grid ρ, and grid V → KS
matrix) — CP2K's "collocation", (3) the facade/`Structure` carrying the `UnitCell` through SCF (§1.D.2), (4)
deciding how GPW's grid-sourced Hartree relates to the existing PW density-fit path (they should share the
FFT-Poisson core, differing only in where ρ comes from). This last point connects to the **future denser-{G} fit
grid** already parked in Round 1 (§6.4/§7 there): GPW's density cutoff is the natural place that lands.

---

## 3. Recommended sequencing

1. **Cheap harmonizations, opportunistically now or as a warm-up:** (C) `IonIon<T>`, and (B) folded into it.
   Self-contained, low-risk, shrink the divergence surface.
2. **`PlaneWave_IBS` bad-habits review (§2.3).** Short, and it de-risks everything after.
3. **Symmorphic space groups (§2.1)**, then **BZ reduction (§2.2a)** as the first, highest-value, GPW-independent
   payoff (validates §2.1 and directly speeds solid SCF). **SALC-PW (§2.2b)** can follow or run in parallel as a
   second validation — not on the GPW critical path.
4. **Facade-preserves-`UnitCell` (§1.D.2)** — small, do it just before / as the first step of GPW.
5. **GPW (§2.4).** The body of the work; it collapses (A) and (D) and is the real target for solids/batteries.

(A) is intentionally *not* a standalone task — it dissolves under GPW. Don't write a molecular
`Integrals_Pseudo<double>`.

---

## 4. Invariants / pins to preserve (carry these into the work)

- **Never assume `orbital == fit`.** Any fit/aux basis is obtained from the orbital basis via
  `Create{CD,Vxc}FitBasisSet(...)` — the factory is the seam even when the answer is trivial (Round 1 §3.1).
- **Fit quality is measured by grid-convergence of ρ, NEVER by ΔE_total** (the fit is non-variational).
- **Spin-polarized is the native formulation**; unpolarized is the ζ=0 collapse. New GPW/periodic terms are
  spin-native from the start (`FittedVxcPol`/`FittedVcorrPol`), unpolarized as the efficiency corner.
- **The principled divergence stays:** real-space `V(r)` (molecular/GPW real-space terms) vs reciprocal `ṽ(G²)`
  (pure PW) is physics, not a wart. Pure plane waves keep `Integrals_Pseudo<dcmplx>` + their intrinsic G-grid.
- **Use well-conditioned bases for any SCF** (Slater/High for atoms; a cleanly-converted GTH valence basis for
  molecular PP) — the "LASolver" symptom is basis conditioning, not a solver gap. `N3`/`N5` are test-only pools,
  invalid for SCF.
- **Regression style:** periodic/GPW energies are "did-E-move" anchors (pin the converged value, no `Converged()`
  guard). Where a real-space-on-lattice quantity must equal its finite counterpart, assert bit-identity
  (`L_PP`-style) rather than an absolute oracle.

---

## 5. Companion documents (history / detail — this doc does not depend on them)

- `doc/MolecularPseudopotentialPlan.md` — the user's PP plan (owned by the user; not edited).
- `doc/MolecularPP_HarmonizationFindings.md` — Round 1: what landed and why (the detailed record §6, the DONE
  divergence work, the fit-basis factory seam, the grid-cutoff analysis).
- `doc/diagrams/pp_molecular_vs_pw.svg` — the map embedded above.
- `doc/FittingCleanupPlan.md` — the fitting-campaign record that delivered Round 1 Item 2 (PW fit-through-factory).
