# `doc/` — what is live, what is a record, what is retired

**Written 2026-09-08** because the answer had stopped being obvious (user: *"Too many plan files, we need
retire some of the old ones"*).  Forty-plus markdown files is not by itself the problem — the problem was
that nothing said which ones a session should READ, which ones it should only CITE, and which ones are
finished.  This file says.  **Keep it current: a doc that changes tier gets its row moved the same day.**

Three tiers, and the rule for each:

| tier | what it means | what a session does with it |
|---|---|---|
| **LIVE** | a QUEUE (`OpenWork.md`, `CleanupCandidates.md`), an INSTRUMENT (`Benchmark.md`, CP2K), or a plan under active execution | READ at session start if you are working in its area |
| **RECORD** | executed; the standing design / clean-code decisions WITH the evidence and the rejected alternatives that produced them | CITE it; never treat its contents as a queue |
| **RETIRED** | `doc/OldPlans/` — fully executed and superseded | historical only |

**Three files, three questions (user ruling 2026-09-16).**  `CLAUDE.md` answers *how do I work here*
(conventions, build/test/box discipline, tool paths, this doc system).  `Pins.md` answers *what must the
code obey* (physics, numerics AND design invariants, one paragraph each, earned by a wrong number).  A
RECORD answers *why is it this way*.  ⇒ **A RECORD holds NO open work.**  When a plan goes RECORD its
residual backlog moves to `OpenWork.md` (or `CleanupCandidates.md` if it is debt), its durable rulings to
`Pins.md`, its conventions to `CLAUDE.md` — the same day.  A "remaining/future" section left inside a RECORD
is a defect of this index; if you find one, harvest it, do not act on it.

---

## LIVE — read these

| file | what it is | state |
|---|---|---|
| **`OpenWork.md`** | ★ **THE tracker.  READ IT AT SESSION START.**  It opens with **THE QUEUED PROGRAMME** — a numbered running order agreed 2026-09-09; start at the first UNFINISHED step (steps 1–4 closed 2026-09-09 → 2026-09-16; **step 5 = DFT+U is next**).  The table under it is the index; everything below that is the evidence that produced it | live |
| **`Pins.md`** | ★ **the durable invariants** — 17 of them: no cut in r space, everything is a fit, integrated observables, spin-native, fit quality = grid-convergence of ρ, … the BasisSet taxonomy (14), smearing/GDM (15), span-matching (16), contemporaneous reporting (17).  **Rulings, not preferences**: each one is there because violating it produced a wrong number at least once | live — read once, then obey |
| **`CleanupCandidates.md`** | the SOLID/OOD debt worklist, and the home of the durable design rulings (R1.0) | live — 34 open items after the 2026-09-08 harvest |
| **`Benchmark.md`** | ★ the standing head-to-head instrument vs CP2K.  **COPY the run command out of §5a; never reconstruct it** | live instrument |
| **`SCFStrategyPlan.md`** | the convergence-acceleration abstraction boundaries (DIIS/GDM, mixing, occupation, loop mode) | ⚠ under the 2026-09-16 rule this is a RECORD candidate (a design note whose one open increment, OT, is already `OpenWork.md` row OT) — not moved, not named in the ruling |
| **`ModuleToolchainPlan.md`** | `import std;` + a modular Blaze fork — banish the preprocessor | ⚠ RECORD candidate under the 2026-09-16 rule (specced, not under execution; a one-line `OpenWork.md` row would carry it) — not moved, not named in the ruling |
| **`LatticeGasPlan.md`** | Li/Na configuration enumeration for the battery work | SPECCED, NOT BUILT — deliberately deferred; ⚠ same RECORD-candidate status as `ModuleToolchainPlan.md` |
| **`BatteryMaterialsRoadmap.md`** | the north star (Li/Na cathode voltage curves) above the individual plans | live orientation |
| **`CP2KBuild.md`** / **`CP2Kresults.md`** | how the oracle is built, and what it says | live reference |

## RECORD — cite, do not queue

| file | what it records |
|---|---|
| **`ParallelAndOraclePlan.md`** | the sequenced road to +U (cut 2026-09-06): Phase 1 our OMP gap ✅ 4.44× (and the lesson — every gain was DELETED serial work, not threading — now in `CLAUDE.md`); 2.1 ✅; 2.5 = programme step 2 ✅; Phase 3 = programme step 5 (its §3.1/3.2 folded into that step); residuals 2.2 + Phase 4 (QE first, on trigger) → `OpenWork.md` row PAR; the OT note → row OT | RECORD — moved 2026-09-16 (user ruling) |
| **`GPWPlan1.md`** | the 2026-07-23 GPW forward queue and the four runtime measurement rounds (2026-08-15→19): the evidence behind the box walk, the stream fold, the BLAS pin.  Its whole pending list was harvested, tree-checked, into `OpenWork.md` Step 3 ("harvested from GPWPlan1") — several items found overtaken; the smearing/GDM finding is pin 15 | RECORD — moved 2026-09-16 (user ruling) |
| **`SphericalLatticePlan.md`** | the 2026-08 MnO accuracy campaign: spherical-d lattice view (I1), the ordering HEALED under it (I2 arm 1) ⇒ **pin 16** (a span can reverse a magnetic ordering); the design rulings (peer implementations behind the abstract capability; the engine is its own dispatch layer).  Its absolute-comparison remainder IS `OpenWork.md` Step 5.  ⚠ its moment conclusions are point-probe numbers (pin 4) | RECORD — moved 2026-09-16 (user ruling) |
| **`BasisSetTaxonomyPlan.md`** | V1.33: **libraries follow the FAMILY (engine), module names carry the GROUP**; §1 is the ruling (now **pin 14**), §3 the library map, §4 the execution log (11 commits, 2026-09-13).  §5's residuals filed: the frozen GPW pairing (OpenWork parked threads), the `qcSymmetry` group-name rename (R2.24), V1.38 | RECORD — row moved 2026-09-16 |
| **`TestSuitePlan.md`** | item **TE** executed in full 2026-09-15: the SCF suite as a checked PRODUCT-SPACE grid — axes (§1), claim kinds (§2), grammar (§3, now a `CLAUDE.md` Tests section), file breakdown (§4), the harness collapse onto `SolidCalculation` + `qchem.Materials` (§6), the `DISABLED_` verdicts (§8, rule now in `CLAUDE.md`), the seven rulings (§11).  Anchor rule → pin 10 | RECORD — row moved 2026-09-16 |
| **`RunReportPlan.md`** | the reporting layer's design (the report IS json; one renderer, layout inferred from shape; a global sink; the section cursor) — migration complete 2026-07-26.  The user's contemporaneous-emission ruling is **pin 17**; the GUI-facing wishlist is ONE `OpenWork.md` parked thread | RECORD — moved 2026-09-16 (user ruling) |
| **`FacadeDFTPlan.md`** | `qchem::Calculation` runs DFT — D1 (unified `Model` enum + `Factory` resolver) and D2 (polarized, via SpinNativeDFTPlan B1–B4) both landed 2026-06-30.  Its "remaining" trio (PBE/GGA, LibXC-polarized, +U) are library increments with rows in `OpenWork.md` (parked threads / programme step 5) | RECORD — moved 2026-09-16 (programme step 4) |
| **`SpinNativeDFTPlan.md`** | spin-native LDA (VWN5, `Ham_DFTcorr_P`, `Molecule_EC(nUp,nDown)`, facade multiplicity) — B1–B4 landed 2026-06-30; the TENET is `Pins.md`, and V1.37 finished the thought (Pol/UnPol = `SpinGroup`) | RECORD — moved 2026-09-16 |
| **`ERI4Rework.md`** | the bra-ket 2× banking — §5 stages 1–3c ALL landed 2026-07-02 (the index's "only 3c remains" was stale for two months; §8 of the file already said DONE).  Its untouched §6 (atomic Rk `LMax`→`Irrep`) is now `CleanupCandidates.md` **R2.23** | RECORD — moved 2026-09-16 |
| **`SCFSeedingPlan.md`** | SAD / IonicSAD / spin-SAD, all built; the "Mn table regen after the d-PP fix" the index carried was DONE 2026-08-06 (`e849b70d`, same day as the fix) — verified against the JSON diff | RECORD — moved 2026-09-16 |
| **`FittingCleanupPlan.md`** | the fitting-layer cleanups; its own header has said *"This file is now a RECORD"* since 2026-09-08 (item C closed as R1.0i) — the index row lagged.  Residuals K (fit-{G} densification, deferred by ruling) and I.1 (the ¾-virial `GetEpsXc` default) live in `CleanupCandidates.md` | RECORD — row moved 2026-09-16 |
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

**Retired 2026-09-14:** `Step2Remaining.md` — the reading aid for the queue's step 2 (rows grouped by what
BLOCKED them).  Every group closed (A/C/D) or was never work (E); the two group-B rows left (**V1.34**,
**R1.0b**) live in `CleanupCandidates.md`.  Worth re-reading for two findings: group D's *"the verdict was
already written in the row"*, and group C's "anchor-moving" rows (V2.2/V2.5) moving nothing when finally
swept.  Its appendix traces the whole `MatrixIntegrator` arc.

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
- **`TestFacadeMigrationPlan.md`** — retired 2026-09-15: both halves executed (molecular 2026-08, solid = TestSuitePlan phase 2)
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
