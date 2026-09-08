# `doc/` — what is live, what is a record, what is retired

**Written 2026-09-08** because the answer had stopped being obvious (user: *"Too many plan files, we need
retire some of the old ones"*).  Forty-plus markdown files is not by itself the problem — the problem was
that nothing said which ones a session should READ, which ones it should only CITE, and which ones are
finished.  This file says.  **Keep it current: a doc that changes tier gets its row moved the same day.**

Three tiers, and the rule for each:

| tier | what it means | what a session does with it |
|---|---|---|
| **LIVE** | has open work in it, and someone is expected to act on it | READ at session start if you are working in its area |
| **RECORD** | executed; kept for the evidence, the durable pins, or a small residual backlog | CITE it; do not treat its contents as a queue |
| **RETIRED** | `doc/OldPlans/` — fully executed and superseded | historical only |

⚠ **A RECORD is not a to-do list.**  Several of these carry a "remaining/future" section that has since been
overtaken; if you find yourself about to act on one, check the LIVE trackers first — they are the authority
on what is actually next.

---

## LIVE — read these

| file | what it is | state |
|---|---|---|
| **`OpenWork.md`** | ★ **THE tracker.  READ IT AT SESSION START.**  The `▶ WHAT IS OPEN` table at the top is the index; everything below is the evidence that produced it | live |
| **`Pins.md`** | ★ **the durable invariants** — no cut in r space, everything is a fit, integrated observables, spin-native, fit quality = grid-convergence of ρ.  **Rulings, not preferences**: each one is there because violating it produced a wrong number at least once | live — read once, then obey |
| **`ParallelAndOraclePlan.md`** | the sequenced phases: Phase 1 (our OMP gap) ✅ closed at 4.44× → Phase 2 (size) → **Phase 2.5 (the SOLID cleanup campaign)** → Phase 3 (DFT+U) → Phase 4 (a second code) | live — this is the road to +U |
| **`CleanupCandidates.md`** | the SOLID/OOD debt worklist, and the home of the durable design rulings (R1.0) | live — 34 open items after the 2026-09-08 harvest |
| **`Benchmark.md`** | ★ the standing head-to-head instrument vs CP2K.  **COPY the run command out of §5a; never reconstruct it** | live instrument |
| **`GPWPlan1.md`** | the GPW forward queue (structs → display → LRU → diffuse+smearing → B_ij) | live |
| **`SCFStrategyPlan.md`** | the convergence-acceleration abstraction boundaries (DIIS/GDM, mixing, occupation, loop mode) | live design note; OT is the open increment |
| **`SphericalLatticePlan.md`** | the MnO accuracy campaign — spherical lattice view, the d-channel, the ordering question | live |
| **`ModuleToolchainPlan.md`** | `import std;` + a modular Blaze fork — banish the preprocessor | live, not started |
| **`LatticeGasPlan.md`** | Li/Na configuration enumeration for the battery work | SPECCED, NOT BUILT — deliberately deferred; the file exists so the design is not re-derived |
| **`BatteryMaterialsRoadmap.md`** | the north star (Li/Na cathode voltage curves) above the individual plans | live orientation |
| **`TestFacadeMigrationPlan.md`** | migrate the molecular tests onto `qchem::Calculation` | live, blocked on the facade gap.  ⚠ See `OpenWork.md` item **TE** first — the test-suite ORGANIZATION question (the axis product, the file breakdown, the naming convention) is the bigger one |
| **`FacadeDFTPlan.md`** | `qchem::Calculation` runs DFT | D1+D2 done; PBE/GGA, LibXC-polarized and +U remain |
| **`ERI4Rework.md`** | the bra-ket 2× banking | stages 1–3b committed; only the 3c cache key remains |
| **`FittingCleanupPlan.md`** | the fitting-layer cleanups | all done except C (the `dynamic_cast` survey) — and that one is Phase 2.5 work |
| **`SCFSeedingPlan.md`** | SAD / IonicSAD / spin-SAD | done through §10; the Mn table entries still need regenerating after the d-PP fix |
| **`SpinNativeDFTPlan.md`** | spin-native XC as the primary formulation | the tenet is live and governs new code (GGA, +U) |
| **`RunReportPlan.md`** | the run report | MIGRATION COMPLETE — but it keeps a real "Remaining / future work" backlog (meta section, detail levels, the `Renderer` DIP split, rolling log) |
| **`CP2KBuild.md`** / **`CP2Kresults.md`** | how the oracle is built, and what it says | live reference |

## RECORD — cite, do not queue

| file | what it records |
|---|---|
| **`GPWPlan.md`** | the 2026-07 GPW campaign record.  ⚠ **Its durable-pins section MOVED to `doc/Pins.md` on 2026-09-08** — nothing in this file is a pin any more.  Its TODO section is superseded by `GPWPlan1.md`.  ▶ Now retirable once its narrative is judged spent |
| **`SymmetryUpgradePlan.md`** | §§0–8 executed; T1/T2/T3 landed and armed; the supercell arc closed 2026-09-08.  Its §9 is a list of DESIGN QUESTIONS to answer when the matching capability is scoped — **not a backlog**.  See its own `▶ STATUS, 2026-09-08` header |
| **`RealComplexPlan.md`** | the real/complex type refactor — the flip is live |
| **`CollocationRewritePlan.md`** | the collocation rewrite — COMPLETE (cache deleted, contract kernel, spin-native XC pair route) |
| **`MolecularPseudopotentialPlan.md`** | Atom_PP + Molecule_PP; the plan that drove the PP interfaces |
| **`MolecularPP_HarmonizationFindings.md`** / **`…Round2.md`** | the molecular↔PW harmonization record.  Round2 §1's "four remaining incidental divergences" is the one live thread in either; ⚠ its own header says GPW has outgrown the file |
| **`GPWGrids.md`** | the table of every grid usage and how its range/spacing is decided |
| **`OTNotes.md`** | what the 2026-07 GDM investigation established, so OT does not re-derive it |
| **`cmakenotes.md`** | build-system notes |

### Histories (append-only closed record; nothing is ever trimmed)

`OpenWork_History1.md`, `OpenWork_History2.md`, **`OpenWork_History3.md`** (cut 2026-09-08),
`CleanupHistory.md`, `GPWHistory.md`, `SymmetryUpgradeHistory.md`, `BenchmarkHistory.md`.

★ **Why they are kept in full**, and it has been earned repeatedly: *a record of what was TRIED AND REJECTED
is worth more than a record of what landed, and it is exactly what gets lost first when a doc is trimmed for
length.*  Three open items in `OpenWork.md` exist only because a measurement refuted the obvious answer.

## RETIRED — `doc/OldPlans/`

Fully executed and superseded.

**Retired 2026-09-08, first pass:** `PlaneWavePlan.md` + `PlaneWavePlan-2.md` (PW-DFT shipped),
`AO_FT_ProjectionCleanup.md` (self-declared DONE, all three moves), `ScreeningPlan.md` (self-declared
*CLOSED — NOT A LIVE PLAN*).

**Retired 2026-09-08, second pass — user corrections to this file's first draft**, each verified against
the tree before moving:

- **`SpaceGroupPlan.md`** — the plan said *workspace `~/Code/qchem7`, branch `lattice-3d-spacegroup`*.
  ⛔ **There is no `qchem7` tree**, and `origin/lattice-3d-spacegroup` is **691 commits behind main**, last
  touched 2026-07-10.  The work LANDED ON MAIN by another route: `src/Symmetry/Lattice_3D/SpaceGroup.C` +
  `Fold.C`, gated by `L_SpaceGroup` / `L_Fold`, and consumed by `UnitCell` / `Lattice_3D` / `BasisSet`.
- **`SymmetryRefactorPlan.md`** — said *IN PROGRESS on branch `symmetry-refactor`*.  ⛔ **That branch does
  not exist**, locally or on the remote.  The `qcSymmetry` reorg landed; `src/Symmetry/{Lattice_3D,Molecule}`
  is its result.
- **`SphericalSALCPlan.md`** — shippable, and its one remainder (**test libcint-spherical, S3b**) is now
  carried by `OpenWork.md` item **TE** as a Stage-C action, which is where it will actually be seen.
- **`APIErgonomicsReview.md`** — it is **GUI-project feedback**, not lib-side work: a separate project built
  on the `pybind/` hooks, and `pybind/` is not ours to edit (CLAUDE.md).  Kept for the record.

⚠ **The lesson worth keeping**: three of these four described a workspace or branch that no longer existed.
A plan file that names a tree or a branch **rots silently** — check `ls ~/Code` and `git branch -a` before
believing one, and prefer *"landed on main as X"* to *"in progress on branch Y"* when writing them.

**Earlier:** `GaussianPlaneWavePlan.md`, `MolecularBasisSetPlan.md`,
`MolecularSymmetryPlan.md`, `ScfrunMoleculeUpgrade.md`, `qcMesh1Design.md`, `qcMeshUpgradePlan.md`,
`FunctionFitterISP.md`, `IntegralCacheUpgrade.md`, `DBCache_issues.md`, `SCF_DIIS_SALC_notes.md`,
`CI_pipline.md`.

## Not prose

`diagrams/` (SVG, embedded in the markdown and in Doxygen), `logs/`, `scripts/`, `Algebra/`, `GSData/`,
`lyx/`, `libcint_ref.pdf`, `Hermite2.ods`.
