# `doc/` — the folder IS the tier

The top of `doc/` is the whole live picture; where a file sits says what it is.  Live docs show only
(1) what is left and how to run the instrument, (2) lessons/gotchas.  History is archived in `Records/`.

**Three files, three questions.**  `CLAUDE.md` = *how do I work here*.  `Pins.md` = *what must the code obey*
(rulings earned by a wrong number).  A RECORD = *why is it this way* (evidence, rejected alternatives).
⇒ **A RECORD holds NO open work**: residual backlog goes to `OpenWork.md` (or the cleanup lists if it is
debt), durable rulings to `Pins.md`, conventions to `CLAUDE.md` — the same day — and the file moves DOWN a level.

⚠ **Source comments cite the old flat paths** (`doc/GPWPlan.md`, `doc/SymmetryUpgradePlan.md`, …); they were
not rewritten (a comment-only sweep forces a module rebuild).  Resolve one with `ls doc/*/<File>.md`.  Relative
links INSIDE moved files (`../src/…`) are one level off.

---

## `doc/` — LIVE (read what applies at session start)

| file | what it is |
|---|---|
| **`OpenWork.md`** | ★ **THE tracker (v4, ≤100 lines).  READ IT AT SESSION START.**  §1 NEXT · §2 MAJOR FEATURES · §3 ACCURACY / OPEN DEFECTS · §4 PERFORMANCE summary.  One row = next concrete action + the one record to read |
| **`HubbardUPlan.md`** | ★ **UNDER EXECUTION** — DFT+U: status, oracle table, ▶ START HERE next actions, gotchas.  Full history: `Records/HubbardUHistory.md`.  Retires to `Records/` when the U-functional lands |
| **`LinearResponsePlan.md`** | **UNDER EXECUTION** — HubbardUPlan's A7: one theory-neutral `ResponseKernel` (DFPT-U / CPHF / MP2 Z-vector), stages R0–R4, each with an oracle.  Retires to `Records/` when A7 lands |
| **`OOD-SOLID-Cleanup.md`** + **`CleanCode.md`** | the cleanup worklist (split 2026-10-01): OOD/SOLID debt vs non-SOLID hygiene, open rows only, original R/V/D ids citing `Records/CleanupHistory2.md` / `CleanupHistory3.md`.  `CleanupCandidates.md` is a 6-line redirect stub |
| **`Pins.md`** | ★ **27 durable invariants** (pin 1: no cut in r space).  Rulings, not preferences; cite as `doc/Pins.md pin N` |
| **`Benchmark.md`** | the standing head-to-head instrument vs CP2K: rules, run commands, the one per-iteration table (§5a), open perf levers (§10).  **COPY the run command; never reconstruct it** |
| **`ModuleToolchainPlan.md`** | `import std;` + modular Blaze fork — banish the preprocessor.  Deferred |
| **`LatticeGasPlan.md`** | Li/Na configuration enumeration.  Specced, not built |
| **`BatteryMaterialsRoadmap.md`** | the north star (Li/Na cathode voltage curves) |
| `README.md` | this index — keep current; a file that changes tier moves folder the same day |

## `doc/Records/` — RECORDS cited by open work (cite; never a queue; nothing here is trimmed)

| file | why it is here |
|---|---|
| **`OpenWork_History1…5.md`**, **`CleanupHistory.md`/`2`/`3`**, **`BenchmarkHistory.md`** (incl. Benchmark.md verbatim as of 2026-10-01), **`GPWHistory.md`**, **`SymmetryUpgradeHistory.md`**, **`HubbardUHistory.md`** | the append-only closed record; trackers WRITE to these (History5 = the v3 tracker verbatim; CleanupHistory3 = the v2 worklist verbatim).  ★ A record of what was TRIED AND REJECTED is worth more than a record of what landed |
| **`CP2KBuild.md`** / **`CP2Kresults.md`** | how the primary oracle is built, and what it says |
| **`SCFStrategyPlan.md`** / **`OTNotes.md`** | convergence-acceleration boundaries; the 2026-07 GDM investigation (prior for row OT) |
| **`ParallelAndOraclePlan.md`** | the road to +U; residuals → row PAR |
| **`GPWPlan1.md`** / **`GPWGrids.md`** | GPW evidence trail; the table of every grid usage |
| **`SphericalLatticePlan.md`** | 2026-08 MnO accuracy campaign (moment numbers are point probes, pin 4) |
| **`SymmetryUpgradePlan.md`** | executed; §9 is design questions, not a backlog |
| **`BasisSetTaxonomyPlan.md`** | the argument behind pin 14 |
| **`TestSuitePlan.md`** | the SCF suite as a checked product-space grid (grammar now in `CLAUDE.md`) |
| **`RunReportPlan.md`** | the reporting layer's standing design (pin 17) |
| **`cmakenotes.md`** | build-system notes for the module toolchain |

## `doc/OldPlans/` — RETIRED (32 files; executed and superseded; historical only)

⚠ **Two lessons the retirements taught:** (1) a plan that names a WORKSPACE or BRANCH rots silently — check
`ls ~/Code` and `git branch -a` before believing one.  (2) an index row saying *"only X remains"* is a CLAIM
about the tree — `git log -S<symbol>` before reading the argument.

## Not prose

`diagrams/` (SVG), `logs/`, `scripts/`, `Algebra/`, `GSData/`, `lyx/`, `libcint_ref.pdf`, `Hermite2.ods`.
