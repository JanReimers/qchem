# Open Work — single tracker for the in-flight threads

One place to see everything dangling, so it isn't spread across sessions. Each item: status, the
next concrete step, and where its plan/code lives. **Recommended discipline: close these before
opening new threads.** Recommended order at the bottom.

Branch: `main` (local, ahead of origin). Anchor build/test: `cd build/Release && ninja UTMain && ./UnitTests/UTMain` (146 green as of this writing).

---

## A. Spherical SALC  ·  PARTIAL (clean checkpoint)  ·  plan: `doc/SphericalSALCPlan.md`

Extend point-group SALC adaptation from the Cartesian PG basis to the two spherical bases.

| Stage | What | Status |
|---|---|---|
| S1 | `SphericalShellRep` (real-spherical operation rep, `qcSymmetry`) | ✅ committed `194f9971`, 4 tests |
| S2 | generalize `AoShell`/`BuildOperationRep` to dispatch Cartesian vs spherical | ⬜ next |
| S3 | the two convention-matched extractors (PG_Spherical in-house; libcint-spherical foreign) | ⬜ the bulk / bug-prone |
| S4 | dispatch in `PG::SymmetryAdapt` by orbital-IBS type | ⬜ |
| S5 | end-to-end tests: spherical-adapted SCF == un-adapted spherical SCF | ⬜ |

**Next step:** S2 — add a spherical angular descriptor to `AoShell` and let `BuildOperationRep` call
`SphericalShellRep` (built S1) for spherical shells; keep the Cartesian path byte-identical.
**Lifts:** the facade's `{.symmetry=true}` Cartesian-only guard, once S3a (PG_Spherical) lands.

---

## B. Spin-native DFT (was "D2 polarized")  ·  NOT STARTED  ·  plan: `doc/FacadeDFTPlan.md` (D2) + tenet `feedback_spin_polarized_primary`

Reframed per the design tenet: **spin-polarized is the native formulation; unpolarized is the
ζ=0 efficiency collapse** — not "add the polarized special case." Four pieces, all spin-first:
1. spin-native VWN5 correlation (ζ-dependent) — *replace* the paramagnetic-only `VWN_Correlation`.
2. `Ham_DFTcorr` carrying both channels (U = collapse) — today only `Ham_DFTcorr_U` exists.
3. open-shell molecular occupation `(n↑,n↓)` — `Molecule_EC(Ne)` is closed-shell aufbau only.
4. facade multiplicity on `CalcOptions`.

**Next step:** write a short spin-native DFT plan doc before code (scope the four pieces). Best
sequenced toward the battery/magnetism payoff, ideally with PBE+U.
**Done already (D1):** closed-shell molecular HF+DFT via the facade — unified `Model` enum + Factory
resolver, LDA + Xα, `8b8df1d0`.

---

## C. Namespace unification (review item 4)  ·  NOT STARTED  ·  plan: `doc/APIErgonomicsReview.md` §4

Move `ScalarFunction`/`Spin`/`Vector3D`/`BasisSet::` under `qchem::`. Whole-tree mechanical sweep;
compile-time-checked (no runtime risk), but huge blast radius.

**Next step:** none yet — **do this LAST**, as one deliberate sweep when the tree is otherwise quiet,
so it never collides with in-flight feature diffs and sweeps the new code in the same pass. Its
consumer benefit is already largely delivered by the `import qchem;` umbrella, so there's no urgency.

---

## D. Test → Facade migration + slim QchemTester  ·  NOT STARTED  ·  plan: `doc/TestFacadeMigrationPlan.md`

Make `qchem::Calculation` THE molecular recipe: migrate the `TestMolecule`-derived molecular tests
(`M_HF_U`, `M_DFT`, `A_DFT` molecular, `M_Sym`→`{.symmetry=true}`) onto the facade, then delete
`TestMolecule` and slim `QchemTester` to the atom/Dirac harness. Atom/Dirac/PP tests stay (facade is
molecule-only). Revives the previously-descoped "fold QchemTester onto the facade" as a deliberate project.

**Prerequisite (additive facade gap):** `CalcOptions` needs `engine`/`angular` selection (libcint /
spherical) for the basis-variant tests — the facade only builds default Cartesian today. (spherical+symmetry
stays guarded until A lands.) Plus extract a free `RelativeError(E,Eref)` test helper. Anchors stay
byte-for-byte; pure refactor.
**Next step:** fill the engine/angular `CalcOptions` gap.
**Sequence:** right after the namespace sweep (so migrated tests use the final namespaces).

---

## Bigger milestone on the horizon (not dangling, just flagged)

- **PBE / GGA** — the highest-value functional for the battery north-star, but a real library
  increment (density-gradient machinery on the mesh), not an enum value. The unified `Model` enum is
  ready to list it with a "not wired" throw. See `doc/FacadeDFTPlan.md`.

## Deferred / descoped (recorded so they're not re-litigated)

- Fold `QchemTester` + the pybind bridge onto the facade (test-harness/binding cleanup; the facade
  makes most of `QchemTester` redundant — a good litmus test, but not lib-surface work).
- Container utils to `src/` (`sample_scalar`/`sample_gradient`, `Structure::BoundingBox`) — binding convenience.
- `SCFParams` ASCII rename — DROPPED; solved instead by C++20 designated initializers (`34ccf302`),
  keeping the compact unicode names.
- `MolecularSym_EC` → `FixedIrrepOcc_EC` rename — belongs with the queued symmetry-naming cleanup.

## Recently closed (for context)

- API ergonomics additive layer: facade (`a7ef787b`), umbrella `import qchem;` (`53206930`),
  `{.symmetry=true}` Cartesian-guarded (`ff5d0a39`), DFT D1 (`8b8df1d0`).
- HF→DFT cross-test energy contamination (fit-basis `Norm` cache not mesh-keyed) — FIXED `c570b20b`.

---

## Recommended order (user-chosen 2026-06-30)

1. **C (namespace sweep)** — NEXT. Tree is quiet now (all committed, A at a clean S1 checkpoint, B/D not
   started), which is exactly when a whole-tree sweep is safe; and doing it before resuming A/D means
   their remaining code is written in the final namespaces. Closes the API-ergonomics core for a GUI-team
   sign-off ("thumbs up when ready", not unblock-ASAP — the umbrella already unblocked them).
2. **D (test → facade migration + slim QchemTester)** — right after C.
3. **A (finish Spherical SALC S2–S5)** — parked at a clean S1 checkpoint; resume when convenient.
4. **B (spin-native DFT)** — plan first, build spin-first, sequence toward PBE+U / batteries.

Note: one library session at a time, GUI on its own branch → land each session at a clean commit and the
whole-tree sweep (C) has nothing to collide with. Hold the line on not opening new threads.
