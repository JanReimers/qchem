# Space-Group Support Plan (`src/Symmetry/Lattice_3D`)

Extend `qcSymmetry` from molecular point groups + Bloch k-labels to crystal
**space groups**. Workspace: `~/Code/qchem7`, branch `lattice-3d-spacegroup`.

Reuses the molecular point-group core (`src/Symmetry/Molecule/`): a space group
≈ point group ⋉ lattice translations, so `SymOp` (the 3×3 orthogonal op),
`AbelianGroup`/`CharacterTable`, `ShellRep`/`OperationRep`, and `BuildSALCs` are
all leverage.

---

## The three sub-tasks (from the spec)

1. **Detect** the ops from cell + basis (spglib-style).
2. **Represent** `{R|τ}` and apply them — to **k** in reciprocal space, to
   **positions/density** in real space.
3. **BZ machinery** — star of k, irreducible wedge (IBZ), integration weights.

## Two tiers

- **Tier A — point group acting on k** (+ time reversal `k→−k`). Consumes only
  the *linear* part `R` of each `{R|τ}`; the fractional translation `τ` is
  irrelevant to how k transforms. Enough for **BZ reduction**. Symmorphic-friendly,
  and correct even for non-symmorphic groups because k reduction never sees `τ`.
  **This is the provider for the qcLattice_3D "BZ-symmetry reduction" item and the
  immediate goal for FCC Si + Monkhorst-Pack.**
- **Tier B — full `{R|τ}` space group.** Fractional translations → non-symmorphic
  `e^{−iG·τ}` phases, density/Hamiltonian symmetrization, crystal SALCs (solid
  analogue of the molecular SALC work), and Wyckoff-position enumeration for the
  GUI. Deferred; the Tier-A data structures are designed not to foreclose it.

---

## How the IBZ interacts with Monkhorst-Pack (the user's question)

The MP grid (`KMesh`) is a uniform sampling: `N = Nx·Ny·Nz` points
`k_ijk = ((i,j,k)+shift)/N` in fractional reciprocal coords, each with weight
`1/N`, approximating `⟨f⟩_BZ = (1/N) Σ_k f(k)`.

Symmetry reduction is **exact, not an approximation**: every k-resolved integrand
we sum (band energies, the k-density contribution `ρ_k`) obeys `f(Rk) = f(k)` for
`R` in the crystal point group `G₀` and `f(−k) = f(k)` by time reversal. So:

1. Partition the `N` grid points into **orbits (stars)** under `G₀` (+ time
   reversal). All points in a star give the same `f`.
2. Keep **one representative** per star. Its weight becomes `|star| / N`.
3. `Σ_reduced w_k f(k) ≡ Σ_full (1/N) f(k)` — identical value, far fewer k-points.

**Grid-closure condition.** Reduction is exact only when `R·k_grid` lands on
another grid point (mod a reciprocal-lattice vector). For a **Γ-centred** grid
this always holds, since `R` permutes the reciprocal lattice and the grid is
`(lattice)/N`. For a **shifted** grid the shift must be symmetry-compatible; if a
particular `R` does not map the grid to itself we simply drop it from the
reduction set (still correct, just less folding) and can warn.

**Index arithmetic.** Point-group ops act on Cartesian k. In the reciprocal
*fractional* basis a lattice symmetry is an **integer matrix** `R̃ = B⁻¹ R B`
(equivalently the transpose-inverse of the direct-lattice integer op). A grid
point `i ∈ [0,N)³` maps to `i' = R̃·i` (for a shifted grid, on `i+shift`), then
reduced mod `N` to find which grid point it is. That integer bookkeeping is the
whole star-finding kernel.

**FCC Si (diamond, Fd-3m #227).** The relevant point group for k-reduction is the
holohedry `m-3m = O_h` (48 ops); centrosymmetric, so time reversal adds nothing.
A 4×4×4 Γ-centred MP grid (64 pts) folds to ~8–10 irreducible points. Note the
*space group* Fd-3m is **non-symmorphic** — but that only matters for Wyckoff/
SALC (Tier B), **not** for the IBZ (Tier A).

---

## Module / folder layout (mirrors `Molecule/`)

```
src/Symmetry/Lattice_3D/
  BlochQN.C                (exists) qchem.Symmetry.Lattice_3D.BlochQN
  SpaceGroup.C             qchem.Symmetry.Lattice_3D.SpaceGroup   -- {R|τ} ops container + detection
  BZReduction.C            qchem.Symmetry.Lattice_3D.BZReduction  -- star of k, IBZ, weights
  (Tier B) CrystalSALC.C   qchem.Symmetry.Lattice_3D.CrystalSALC
  Imp/
    SpaceGroup.C
    BZReduction.C
```
Namespace `qchem::Symmetry::Lattice_3D`. Registered in
`src/Symmetry/CMakeLists.txt` alongside the existing `Lattice_3D/BlochQN.C`
(PRIVATE `Imp/*.C`, PUBLIC module-interface `*.C`).

---

## Incremental delivery

### Increment 1 — Tier-A IBZ reduction for FCC Si (the "first" deliverable) — DONE (uncommitted)
Pure symmetry math, **not gated on the complex-D SCF bug** (below). Unit-testable
in isolation. Decided with user: scope = IBZ math + tests only (no GPW wiring yet);
general holohedry detector (not FCC-focused).

Delivered:
- `qchem.Symmetry.Lattice_3D.SpaceGroup` — `Detect(A, basis, tol)`; holohedry (integer
  ops `WᵀMW=M`) ∩ basis symmetry (∃τ); `PointGroupOps`, `ReciprocalPointOps`, `isSymmorphic`.
  Takes cell-matrix + fractional `AtomSite`s (qcStructure-free, mirrors `Molecule/PointGroup`).
- `qchem.Symmetry.Lattice_3D.BZReduction` — `ReduceToIBZ(N, shift, reciprocalOps)` →
  `IBZMesh{points, ownerOfGrid}`; orbit-closure BFS; lowest-index representative per star.
- Tests `src/Symmetry/tests/L_SpaceGroup.C` (58/58 UTSymmetry green): FCC Si → O_h 48 ops
  non-symmorphic; SC 2×2×2 Γ → 4 reps {1,3,3,1}; **FCC 4×4×4 Γ → 8 irreducible k-points**;
  identity-ops no-fold. `UTMain` link-builds (no downstream break).

Original sub-steps:

1. `SpaceGroup` (linear part only for now): detect the crystal point group `G₀`
   from `UnitCell` = holohedry of the Bravais lattice (integer ops `M` with
   `Mᵀ G M = G`, `det = ±1`, `G` = metric) ∩ ops that map the atomic basis to
   itself (allowing an existence-of-`τ` check). Reuse `SymOp` for the Cartesian
   3×3 rep; hold the integer lattice rep `R̃` for the grid arithmetic.
2. `BZReduction`: given a `KMesh` (unreduced MP) + `G₀`, compute stars, pick
   representatives, aggregate weights → reduced `{k}`, `{w}`.
3. Expose `Lattice_3D::MakeIBZMesh(shift)` (or an overload) returning the reduced
   `KMesh`, so the GPW/PW k-loop is unchanged — weights already flow through
   `BlochFactory(N,ik,weight,shift) → BlochQN::GetWeight() → TOrbitalsImp::GetChargeDensity`.
4. **Tests** (combinatorial, no SCF): weight sum = 1; FCC Si O_h on 2×2×2 / 4×4×4
   Γ-centred grids gives the textbook irreducible count + weights; `Σ_reduced w·k-orbit`
   reproduces the full grid; a low-symmetry cell reduces trivially.

### Increment 2 — wire IBZ into GPW and validate the energy
Swap `MakeKMesh` → `MakeIBZMesh` in the GPW k-loop; assert reduced-mesh energy
== full-mesh energy for FCC Si. **Gated on GPWPlan TODO 1 (complex-D SCF bug)** —
interior irreducible k-points have genuinely complex Bloch phases, which currently
over-bind and fail to converge. Increment 1 stands alone without this.

### Increment 3 — Tier B: full `{R|τ}`, crystal SALCs, Wyckoff (GUI)
Fractional translations, `e^{−iG·τ}` phases, per-k-block SALCs (reusing
`BuildSALCs`), density symmetrization, and Wyckoff enumeration for GUI
space-group entry. **Wyckoff needs Tier B** (positions are orbits under the full
`{R|τ}`; diamond's site symmetry comes from the non-symmorphic ops).

---

## Notes / risks

- **From-scratch detection**, reusing the `Molecule/` core, consistent with the
  project's build-from-scratch symmetry lineage (no spglib dependency). spglib
  remains a possible oracle for cross-checking in tests.
- **Complex-D bug** (GPWPlan.md TODO 1) gates the *energy validation* (Increment
  2), not the reduction math (Increment 1).
- **GDM/Ladder are `<double>`-only** — affects convergence speed of complex-k SCF
  even after the complex-D fix.
- Carry per-k weights from the start (already true) so reduction folds in without
  restructuring the k-loop.
