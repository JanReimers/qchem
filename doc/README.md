# `doc/` — the folder IS the tier

**Rewritten 2026-09-16** on the user's ruling: *"I would like to keep these in doc folder: OpenWork.md Pins.md
CleanupCandidates.md Benchmark.md ModuleToolchainPlan.md LatticeGasPlan.md BatteryMaterialsRoadmap.md.  It
would help me if most of the others are retired … or if you find any of them useful for the pending
Open/Cleanup/Toolchain you can put those in a Records folder."*  So the top of `doc/` is now the whole live
picture, and where a file sits says what it is.

**Three files, three questions.**  `CLAUDE.md` answers *how do I work here* (conventions, build/test/box
discipline, tool paths, this doc system).  `Pins.md` answers *what must the code obey* (physics, numerics
AND design invariants, one paragraph each, earned by a wrong number).  A RECORD answers *why is it this way*
(the evidence and the rejected alternatives).  ⇒ **A RECORD holds NO open work**: when a plan is executed
its residual backlog moves to `OpenWork.md` (or `CleanupCandidates.md` if it is debt), its durable rulings
to `Pins.md`, its conventions to `CLAUDE.md` — the same day — and the file moves DOWN a level.  A
"remaining/future" section found inside a record is a defect of this index: harvest it, do not act on it.

⚠ **Source comments cite the old flat paths** (`doc/GPWPlan.md`, `doc/SymmetryUpgradePlan.md`, …).  They
were deliberately not rewritten — a comment-only sweep of module interface units forces a rebuild.  Resolve
one with `ls doc/*/<File>.md`.  Relative links INSIDE moved files (`../src/…`) are likewise one level off.

---

## `doc/` — LIVE (eight files; read what applies at session start)

| file | what it is |
|---|---|
| **`OpenWork.md`** | ★ **THE tracker (v3, rebuilt 2026-09-16).  READ IT AT SESSION START.**  §1 NEXT (step 5, DFT+U) · §2 MAJOR FEATURES (DFT+U, GGA, hybrids, forces, USPP, OT, k-parallelism, space-group irreps, SSB descent, …) · §3 NON-OOD CLEANUP · §4 REMAINING TODO (accuracy, performance) · parked · descoped.  One row = next concrete action + the one record to read |
| **`CleanupCandidates.md`** | the SOLID/OOD debt worklist (v2, rebuilt 2026-09-19): the user's charter + ~30 open rows in R/V/D tables, each citing `Records/CleanupHistory2.md` by id; the ~70 closed rows and every ruling's full argument live there |
| **`Pins.md`** | ★ **23 durable invariants** — no cut in r space (1), everything is a fit (2), integrated observables (4), spin-native (5), … the BasisSet taxonomy (14), smearing/GDM (15), span-matching (16), contemporaneous reporting (17), the XC feed / mixer selectivity (18), never D-screen the gather (19), ask what a matrix means (20), pivoted Cholesky + canary (21), vet-stage equivariant trim (22), +U is orbital-resolved and the manifold is an input (23).  Rulings, not preferences; cite as `doc/Pins.md pin N` |
| **`Benchmark.md`** | ★ the standing head-to-head instrument vs CP2K.  **COPY the run command out of §5a; never reconstruct it** |
| **`ModuleToolchainPlan.md`** | `import std;` + a modular Blaze fork — banish the preprocessor.  Deferred, not started |
| **`LatticeGasPlan.md`** | Li/Na configuration enumeration for the battery work.  Specced, not built — kept so the design is not re-derived |
| **`BatteryMaterialsRoadmap.md`** | the north star (Li/Na cathode voltage curves) above every plan |
| `README.md` | this index.  **Keep it current: a file that changes tier moves folder the same day** |

## `doc/Records/` — RECORDS still cited by open work (cite; never a queue)

| file | why it is still here |
|---|---|
| **`OpenWork_History1/2/3/4.md`** (4 = the v2 tracker verbatim, 2026-09-16), **`CleanupHistory.md`** + **`CleanupHistory2.md`** (2 = the v1 worklist verbatim, 2026-09-19), **`BenchmarkHistory.md`**, **`GPWHistory.md`**, **`SymmetryUpgradeHistory.md`** | the append-only closed record; the trackers WRITE to these.  ★ *A record of what was TRIED AND REJECTED is worth more than a record of what landed* — nothing here is ever trimmed |
| **`CP2KBuild.md`** / **`CP2Kresults.md`** | how the primary oracle is built, and what it says; the +U oracle row (step 5) starts here |
| **`SCFStrategyPlan.md`** / **`OTNotes.md`** | the convergence-acceleration abstraction boundaries, and what the 2026-07 GDM investigation established — the design and the prior for row **OT** |
| **`ParallelAndOraclePlan.md`** | the sequenced road to +U: Phase 1 ✅ 4.44×, 2.5 = programme step 2 ✅, Phase 3 = step 5 (folded in), residuals → row **PAR** |
| **`GPWPlan1.md`** | the GPW forward queue's evidence trail (four runtime rounds); its pending list is harvested into `OpenWork.md` Step 3 |
| **`SphericalLatticePlan.md`** | the 2026-08 MnO accuracy campaign; its absolute-comparison remainder IS `OpenWork.md` Step 5; pin 16.  ⚠ its moment conclusions are point-probe numbers (pin 4) |
| **`SymmetryUpgradePlan.md`** | §§0–8 executed (T1–T3, the supercell arc, the MnO Shubnikov campaign); its §9 is DESIGN QUESTIONS for when the matching capability is scoped, not a backlog; rows **KP**, **BM**, R1.0r cite it |
| **`BasisSetTaxonomyPlan.md`** | V1.33: libraries follow the FAMILY, modules carry the GROUP — the argument behind pin 14, and the parent of V1.38 / R2.24 |
| **`GPWGrids.md`** | the table of every grid usage and how its range/spacing is decided — item 1 (size the Becke grid) reads it |
| **`TestSuitePlan.md`** | item TE: the SCF suite as a checked product-space grid; the grammar is now in `CLAUDE.md`, `scripts/testgrid` checks it; remainders (S3b, PW facade axis, `_Long` budget) are tracker rows |
| **`RunReportPlan.md`** | the reporting layer's standing design (the report IS json; one renderer; a global sink; the section cursor) — pin 17's record; the GUI wishlist is a parked thread |
| **`cmakenotes.md`** | build-system notes — `ModuleToolchainPlan.md` will need them |

## `doc/OldPlans/` — RETIRED (executed and superseded; historical only)

Thirty-two files.  **Retired 2026-09-16** on the ruling above: `GPWPlan.md` (the 2026-07 campaign; its pins
left for `Pins.md` on 2026-09-08), `RealComplexPlan.md` (track complete 2026-08-19), `CollocationRewritePlan.md`
(complete 2026-08-28), `SpinNativeDFTPlan.md` (B1–B4 landed 2026-06-30; the tenet is pin 5), `FacadeDFTPlan.md`
(D1+D2 2026-06-30), `SCFSeedingPlan.md` (all phases + spin-SAD; the Mn regen was done 2026-08-06),
`FittingCleanupPlan.md` (residuals K + I.1 are `CleanupCandidates.md` rows), `ERI4Rework.md` (§5 landed
2026-07-02; §6 is R2.23), `MolecularPseudopotentialPlan.md`, `MolecularPP_HarmonizationFindings.md` +
`…Round2.md` (the molecular↔PW harmonization, outgrown by GPW).

**Retired earlier** (2026-09-08/14/15, each verified against the tree first): `Step2Remaining.md`,
`PlaneWavePlan.md` + `-2.md`, `AO_FT_ProjectionCleanup.md`, `ScreeningPlan.md`, `SpaceGroupPlan.md`,
`SymmetryRefactorPlan.md`, `TestFacadeMigrationPlan.md`, `SphericalSALCPlan.md`, `APIErgonomicsReview.md`,
`GaussianPlaneWavePlan.md`, `MolecularBasisSetPlan.md`, `MolecularSymmetryPlan.md`, `ScfrunMoleculeUpgrade.md`,
`qcMesh1Design.md`, `qcMeshUpgradePlan.md`, `FunctionFitterISP.md`, `IntegralCacheUpgrade.md`,
`DBCache_issues.md`, `SCF_DIIS_SALC_notes.md`, `CI_pipline.md`.

⚠ **Two lessons the retirements taught, kept here because they recur:** (1) a plan that names a WORKSPACE
or a BRANCH rots silently — three of the 2026-09-08 retirements described trees or branches that no longer
existed; check `ls ~/Code` and `git branch -a` before believing one.  (2) an index row that says *"only X
remains"* is a CLAIM about the tree — on 2026-09-16 five such rows were checked and **none** was open;
`git log -S<symbol>` before reading the argument.

## Not prose

`diagrams/` (SVG, embedded in the markdown and in Doxygen), `logs/`, `scripts/`, `Algebra/`, `GSData/`,
`lyx/`, `libcint_ref.pdf`, `Hermite2.ods`.
