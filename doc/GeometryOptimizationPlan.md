# Geometry optimization: the molecular start (DESIGN NOTE, step 1)

Lane B, 2026-11-07.  Status: a design for the user to read before any code.  It covers the OpenWork rows "A cheap host for
geometry-optimization experience" (the user's sequence: molecules, then rutile TiO₂, LiCoO₂, Li_xTiS₂) and level 1 of "Geometry
optimization: imposed relaxation, then DISCOVER…".  It builds on `doc/Records/CalcNetworkReview.md` §3 (the
`PotentialEnergySurface` face) and on the seed test (`doc/Records/SeedTestReport.md`).  The measurements here are a PySCF probe
on H₂O, `~/Code/qchem6-runs/H2O/fdforces/` (`fd_h2o.py`, `dh_h2o.py`, their `.json`).

## 0. What exists, and what this note adds

- **Already there.**
  - `qchem.PotentialEnergySurface` holds `ConvergeAt(cell, R, want) → Outcome<PESPoint, PESFailure>`.
  - Two implementers exist. `MoleculeSession` carries each point's D block by block onto the next point. `SolidSession` drives the
    deck's `displace` action and sizes its grids on max|ΔH_ij| (`ClientConsumes::MatrixElements`).
  - A forces request returns `NoGradient` today.
- **Added here.**
  - A finite-difference forces layer.
  - A structure-neutral optimizer that sees only that face.
  - The rule for the grids along an optimization.
  - The gaps in the molecular facade.
  - The error budget that ties the FD step h to ε.
  - What pin 32 means for FD forces now and for analytic forces later.

```
   GeometryOptimizer ──drives──► PotentialEnergySurface ◄── FiniteDifferenceForces (decorator: forces from 6N energies)
     (BFGS + trust radius,          (energy, forces)              │ wraps any
      knows only R, E, F)                                         ▼
                                                     MoleculeSession │ SolidSession   (one fixed setup each)
```

## 1. Finite-difference forces through `PESPoint` energies

`FiniteDifferenceForces : PotentialEnergySurface` wraps another surface.
- A `ConvergeAt(…, {forces=true})` on it converges the centre, then 6N displaced points (central differences, ±h along each
  Cartesian).
- It returns the centre's `PESPoint` with `forces = −(E₊ − E₋)/2h`.
- A request without forces passes straight through.

Why this shape:
- It is structure neutral: it never asks what is under the face.
- It keeps the optimizer ignorant of where forces come from.
- When analytic forces land it stays as their test oracle, at any point on either surface. That is the role §3.3 of the review
  already gives FD.

Ordering. The inner session carries the LAST point's D. So the decorator visits +h then −h per coordinate, which is a 2h step,
and ends at the centre's neighbourhood; the next optimizer step then starts from a D within h of where it was converged.

Each point's iteration count stays in the report, as a `report::Row` per displacement under the step's row.

Symmetry: a Cartesian displacement breaks an imposed group (`SymmetryBroken`), so the molecular start runs with
`symmetry = false`. Displacing along symmetry-adapted coordinates is level 1's business (§7).

**Cost.**
- H₂O: 18 SCFs per gradient. Each is a carried restart of a few iterations (our water SCF at the default mesh is about 3 s
  single-threaded cold).
- So a gradient is under a minute, and a 6-10 step optimization is minutes.
- Benzene (12 atoms, 72 SCFs per gradient) is the practical ceiling for FD.

**Oracle.** PySCF gradients, same basis, PP and functional. On H₂O, PySCF's analytic gradient (level-9 grid, grid response on)
and our central FD on the (250, 59) reference mesh agree to 3.5e-7 Ha/bohr (stretch 0.0153565 vs 0.0153568). That is the gate
resolution available.

Gates, after the code (§8): H₂O, NH₃, CH₄ at a displaced geometry, with FD forces against PySCF gradients to the §5 budget. Then
the optimized geometry against PySCF's own optimizer (geomeTRIC or PyBerny) to 1e-4 bohr, plus our step count.

## 2. The optimizer: Cartesian BFGS with a model Hessian and a trust radius

**Recommendation.**
- Cartesian quasi-Newton (BFGS update of the Hessian, rational-function or trust-region step), started from a model Hessian
  (Lindh's force-field form, which needs only Z and R).
- External modes are projected out: 3 translations and 3 rotations for a molecule (cell empty), 3 translations for a solid.

**Why not internal coordinates.**
- Redundant internals (Pulay/Fogarasi, geomeTRIC's TRIC) cut molecular step counts by about 2×. But they are a molecular
  construction, and they are the hard part of an optimizer.
- The user's sequence leaves molecules after the first exercise, and in a crystal the natural coordinates are the free Wyckoff
  parameters and the cell (level 1), not bonds and angles.
- With the Lindh Hessian, Cartesian BFGS keeps most of the gain and stays structure neutral.
- This is a judgement, not a measurement: the H₂O/NH₃/CH₄ step counts against geomeTRIC will show the gap.

**Shape.**
- `GeometryOptimizer(PotentialEnergySurface&, OptimizerOptions)::Run(start) → Outcome<OptimizedGeometry, OptimizeFailure>`.
- No line search: a trust radius (start 0.3 bohr), shrunk on an energy rise, with that step rejected.
- A `PESFailure` at a step is the optimizer's outcome too (stop, report which point and why), never a throw.

**Convergence** is max|F|, RMS F, max step and ΔE, the usual four. The force threshold is derived from ε, never stated in
isolation (§5): a threshold below the force noise floor makes the optimizer chase the grid.

**Report.** One `report::Row` per step: E, max|F|, step length, trust radius, SCF iterations summed over the step's points.

## 3. Grids fixed across one optimization, picked once from the two-step run

The rule (the user's fixed-settings and one-grid-set-per-series rulings, applied to one optimization):

1. **Step 1.** Seed (SAD), ε = 1e-4, converge at the start geometry.
2. **Audit at that D on the H measure** (max|ΔH_ij| per rung against a fixed fine reference, the persistent pick), at the ε the
   optimization's force target needs (§5).
3. **Fix the RULE for the whole optimization:** (nR, L, angular kind, radial map, Becke order) per species, the Hartree fit basis
   and every raster.  The mesh POINTS move with the atoms (atom-centred); its parameters never change.  The basis recipe and the
   PP set are frozen too, as `MoleculeSession` already freezes them.
4. **Never resize mid-optimization.**  A resize makes E(R) discontinuous and BFGS's Hessian update meaningless.
5. **Report health per point:** the mesh's charge error δq (free at every point).  At the converged geometry, run the same audit
   again and report whether its pick differs.  A different pick is reported against ε, and acting on it is the NEXT run's
   business (pin 31 catch 2: a miss is acted on at a run boundary).

Evidence:
- The seed test: on all five systems the 1e-5 pick made at the 1e-4 D equals the converged D's pick (0 rungs, ×1.00 points),
  and the seed alone is unreliable at 1e-5.
- What that test did NOT cover is the geometry change between the two steps. The final-geometry audit of step 5 is the check
  that closes it, for every optimization rather than once.
- The audit's cost (on MnO the full ladder is 15× the SCF) is paid once per optimization, so on a molecule it is small next to
  the 6N-point gradients.

## 4. What the molecular facade lacks today

| gap | today | needed for the H₂O start |
|---|---|---|
| ε | `CalcOptions::accuracy` exists but is consumed by no term ("consumed by: none yet"); ε is wired on the solid facade only (V-FACADE-ACCURACY) | the molecular XC mesh (and Hartree fit route) sized from ε on H; at minimum the audit below, with the pick fed back as the setup's `mesh` |
| a molecular grid audit | none in our code: the seed test's TeF₆ ladder was PySCF; the solid audit is `SolidCalculation::XCMeshAudit` | `MoleculeCalculation` audit at a fixed D: E_xc and the V_xc matrix per rung against a reference, normalised-basis max\|ΔH_ij\| |
| two steps across setups | a `MoleculeSession` carries D between points of ONE setup | carry D from the loose setup's converged point into the tight setup (same basis, so D carries; only the mesh differs) |
| the deck | refuses molecules ("the deck currently drives solids") | molecules in the deck, and an `optimize` post-action beside `displace` |
| basis | the transcribed DZVP-MOLOPT `.bsd` is H C N O only (gated against CP2K, 2026-11-02) | enough for H₂O/NH₃/CH₄; `scripts/cp2k2bsd` for Li, Ti, Co, S, F when the solids come |
| forces | `NoGradient` | §1's decorator now, analytic later (§6) |

## 5. The FD step h against the grid: what ε must be

Three errors add in an FD force:

    F_FD − F = ∇δE_grid(R)  +  (h²/6) E'''  +  O(τ_SCF / h)

Measured on H₂O (LDA, GTH, DZVP-MOLOPT, our molecular rule MHL m3 α2 Becke order 2, force on one H):

| rung (nR, L) | points | \|ΔE\| at R₀ (μHa) | max\|ΔH_ij\| | grid force error, stretch / bend (Ha/bohr) | force/H |
|---|---|---|---|---|---|
| (40, 17) | 13,164 | 1.2 | 2.0e-4 | 2.1e-3 / 6e-7 | 10 |
| (60, 23) | 34,864 | 2.0 | 6.8e-6 | 9.5e-5 / 3.0e-5 | 14 |
| **(80, 29) default** | 72,282 | 0.09 | 3.6e-7 | 4.9e-6 / 2.7e-6 | 14 |
| (100, 35) | 129,722 | 0.13 | 6.2e-8 | 4.5e-7 / 4e-8 | 7 |
| (150, 47) | 345,292 | 0.004 | 2.3e-9 | ~1e-8 | — |

What it shows:

1. **On an atom-centred mesh that moves with the atoms, the grid error is a smooth BIAS, not noise.**
   - The rung's force error barely moves with h: within 1% from h = 1e-3 to 1e-2 bohr, and 8% at 3e-2 (e.g. (80,29): −4.93,
     −4.92, −4.88, −4.53 μHa/bohr).
   - FD differentiates the grid's own δE(R) exactly, so the grid does not set h.
   - (A uniform raster that does NOT move with the atoms is different: its egg-box error oscillates with the raster spacing, and
     there h and the grid interact. That is the solid case, TiO₂ onward.)
2. **The force error follows max|ΔH_ij|, not ΔE.**
   - (40, 17) has the smaller energy error of the first two rungs and the 20× larger force error.
   - Across the rungs the force error is 7-14× max|ΔH_ij|. That is ONE molecule and two coordinates, an empirical ratio to
     calibrate on more systems, never a constant to fold in (CLAUDE.md, rigorous vs empirical).
   - This is pin 31's measure clause, now measured on a force.
3. **Truncation.**
   - Stretch: (h²/6)E''' gives 0.3, 2.8, 34 and 312 μHa/bohr at h = 1e-3, 3e-3, 1e-2, 3e-2 bohr.
   - The bend is flat to 1e-8.
4. **SCF tolerance.**
   - At h = 1e-3 an SCF converged to ΔE 1e-8 adds 7e-10 Ha/bohr, but one at 1e-6 adds 2.7e-4.
   - The energy must be converged to ≤ 1e-9 Ha relative between the ±h points. Our SCF gates on [F,D] and Δρ, so the FD layer
     must state a tight `MinFD` for its points (the E error is quadratic in the residual for a consistent Hamiltonian, §6).

So, for a force target f_tol (Ha/bohr):
- **h = 2e-3 bohr** (1e-3 Å): truncation ≈ 1 μHa/bohr.
- **SCF to 1e-9 Ha.**
- **ε_H ≈ f_tol / 15** on the audit.

Examples:
- The conventional "default" threshold (max F 4.5e-4 Ha/bohr) needs ε_H ≈ 3e-5, so today's ε = 1e-5 meets it with margin.
- A "tight" 1.5e-5 Ha/bohr needs ε_H ≈ 1e-6, a finer rung than the 1e-5 default delivers.
- The molecular default mesh (80, 29) delivers about 5e-6 Ha/bohr on H₂O.

Proposal:
- The optimizer's default force threshold is derived as 15·ε_H × a margin of 3, so 4.5e-4 at ε = 1e-5. That is the
  conventional default by coincidence, not by design.
- The run prints the rigorous part (none: no bound exists for the XC mesh), the calibration (15, one system) and the result
  side by side, in one place.

## 6. Pin 32: FD now, analytic later

- **FD forces differentiate the E the run reports, whatever H is.**
  - If a term's H were not ∂E/∂D, the SCF would stop where [F, D] = 0 but E is not stationary in D. E(R) would then carry a
    first-order D error, and FD would still return the true derivative of that E(R), the surface the optimizer actually
    samples. So FD forces are right for the surface they sample.
  - The price is noise: E is no longer quadratic in the SCF residual, so §5's SCF tolerance would have to tighten.
  - The default molecular stack is consistent: Hartree = `Vee_DFT_Dunlap` (Cᵀ G⁻¹ C) and Vxc = the δ fit on the mesh's points.
    Both are on pin 32's consistent list.
  - The A1 Vxc fit (`VxcFit::Auxiliary`) and any overlap-metric route must be checked with pin 32's test before an
    optimization uses them.
- **Analytic forces need H = ∂E/∂D for every term.** Hellmann-Feynman + Pulay (the review's §3.3 formula) holds only at a
  stationary point of E in D, so each fitted term must pass the symmetric-second-derivative test first. Then come each term's
  explicit R-derivative pieces:
  - the moving fit basis (Dunlap's ∂C/∂R, ∂G/∂R);
  - the Becke weight derivatives (the mesh moves);
  - the PP's projector derivatives;
  - the Pulay term −Tr(W ∂S/∂R).

  The FD layer of §1 is each piece's oracle.  The analytic route is the forces note's (lane A, C3); this note leaves the face
  unchanged for it.

## 7. Toward the solids (out of scope here, named so the molecular design does not block them)

- **TiO₂ (level 1, imposed):** the variables are a, c and u.
  - A cell change is `CellChange` today, so either stress or FD on the cell is needed. FD on the cell needs a setup whose grids
    follow the cell by RULE (the series ruling's "same rules"), not by count.
  - Wyckoff u by FD needs a symmetry-preserving displacement: the decorator gains "displace along these symmetric modes"
    instead of Cartesians.
- **Uniform rasters** bring the egg-box: §5's "bias, not noise" fails there, and h must be chosen against the raster spacing.
  Measure on TiO₂ before trusting FD forces on a solid.

## 8. Steps and estimate

| step | content | estimate |
|---|---|---|
| 1 | this note (user review) | done |
| 2 | `FiniteDifferenceForces` decorator + `GeometryOptimizer` (BFGS, Lindh Hessian, trust radius, external modes projected); unit tests on an analytic model surface (Morse / LJ cluster through a test PES, no SCF) | 1 day |
| 3 | the molecular audit (E and normalised H per rung at fixed D, reported), the two-step pick, and the D carry across setups | 1.5-2 days (overlaps V-FACADE-ACCURACY and the ε inventory: coordinate with lane S) |
| 4 | gates: H₂O/NH₃/CH₄ FD forces vs PySCF gradients; optimized geometries vs PySCF's optimizer, step counts; ITMain, `_Long` if needed | 1 day + minutes of runs |
| 5 | molecules in the deck + `optimize` post-action | 1-2 days (the deck refuses molecules entirely today) |

About 5-7 lane-days to an H₂O optimization from a deck on FD forces, with grids fixed from the two-step audit.  Rutile TiO₂
follows with the cell and egg-box questions of §7.

## 9. Questions for the user

1. **Coordinates:** Cartesian BFGS with a model Hessian (structure neutral, recommended), or redundant internals for molecules
   (fewer steps, molecular only)?
2. **The force threshold:** derive it from ε (§5's 15·ε_H × 3, one-system calibration), or state it, with ε derived from it?
3. **FD's reach:** is benzene-size (≈ 70 SCFs per gradient) enough for the molecular phase, so analytic forces wait for the
   forces note?
