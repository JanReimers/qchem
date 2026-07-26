# RunReportPlan — a professional, JSON-ready reporting layer

Status: PLAN (2026-07-26).  Consolidate the scattered setup diagnostics into one organized, machine-readable
report — rendered as tidy terminal output today and JSON for the GUI tomorrow, from a single data model.

## Motivation

The information a user needs to understand a run already exists, but as ad-hoc `cout` scattered across
subsystems and printed in whatever order the code happens to reach it:
- `[GPW grid] ladder L=… N=… Ecut=… nG=…` (per collocation level) — `BasisSet/Lattice_3D/Evaluators/GPW`.
- `[overlap S] n=… min eig=… min sv=… max eig=… cond=…` — `LASolver` (`report_conditioning`, toggle-gated).
- `[ortho] Auto: overlap has k near-null mode(s) …` — `LASolver`.
- `IntegralsCache RAM usage report … PG.Omega … reuse=… RAM=…` — the evaluator cache stats.
- SCF settings — currently NOT reported at all (buried in `SCFParams`); the reader guesses the recipe.

None of it lines up in columns, the order is nondeterministic, and none of it is machine-readable.  The GUI
app (nanobind + PySide6 + pyqtgraph + HDF5 — see `project_viz_gui_plan`) will need the SAME data as JSON.

## Core principle: the report IS json

Instead of hard-coded C++ section structs, take the data model all the way to `nlohmann::json` itself, and
separate the report's **data** from its **rendering** (one of the anticipated renderers just happens to also
BE json).  Each **data provider** (basis set, `LASolver`, GPW grids, …) decides what to include, and the
content extends without touching a struct definition.  We give up compile-time type safety — but a report is
**display, not physics**: a wrong key is a cosmetic bug, not a wrong number, so the project's "prefer
compile-time" bias (which rightly governs the wavefunction types) is the wrong tool here.  **Typed physics,
schemaless display.**  And data providers/consumers only need to be *partially* in sync: a provider may always
supply MORE than a consumer reads (additive schema) — the sync is one-directional and loose.

```
report DATA  ==  an nlohmann::json          (the model; leaf module qchem.Reporting -- depends only on json + Types)
     |           data providers WRITE it key-free (into CurrentRunReport; see "Ownership")
     |
     +--> RenderConsole(const json&, ostream&)   -> terminal, today  (ONE generic renderer -- table or tree)
     +--> RenderJson(const json&)                -> report.dump()    (GUI/HDF5, tomorrow; trivial)
```

## Rendering: one generic renderer, layout INFERRED from json structure

The layout is NOT annotated with hints — it is **inferred from the json's shape**.  A provider already
communicates its intent by HOW it shapes the data (an array of uniform objects means "table"; a nested object
means "tree"), so every renderer — console now, GUI later — infers the same layout from the same structure.
No presentation metadata ever pollutes the data json.  (If some pathological section ever truly needed an
explicit hint, it would go in a SEPARATE metadata json, never mixed into the data — but we expect never to.)

Only **two layouts**, chosen by interrogating a section's nesting depth (and breadth):
- **depth ≥ 3** (objects nested within objects) → **tree**.
- **depth 2** (array/object of flat objects) → **table** if rows ≥ 2, else **tree** (a 1-row table reads worse
  than a key/value block).  The row threshold is a tunable knob — start at "≤1 → tree".
- **depth 1** (object of scalars — `scf.standard`, `meta`) → the degenerate **two-column key/value table**.

Tabulate renders the tables; a simple indented walker renders the trees.  **We push HARD to keep the renderer
generic — NO per-section renderers.**  Sections that don't tabulate (the eigen spectrum, the cache summary)
are simply SHAPED as trees, and the generic tree path prints them; we accept console-layout compromises to
keep the renderer generic.  A per-section override survives only as a theoretical last resort for a section
that is neither table nor tree (e.g. a console sparkline/plot) — and the goal is to never need it.

`RenderConsole` is non-trivial (it may be a small class holding the generic table + generic tree walkers), but
it stays GENERIC.  `RenderJson` is `report.dump()`.

## Detail level (verbosity) — a CONSOLE-ONLY concern

Like every logging framework, console output has a **detail level** (Terse / Normal / Verbose).  The rule:
**data providers ALWAYS store everything, regardless of the level** — the report (hence the JSON / GUI / HDF5
record) is ALWAYS complete.  The level is consumed ONLY by `RenderConsole`, which decides how much of the
complete report to SHOW.  Log everything; filter at output.  The "what shows at which level" map lives on the
render side (a section/optional-field minimum level), so the json never becomes a half-populated record just
because someone ran Terse.  Note this is a DIFFERENT axis from layout: layout (table/tree) comes from
structure; verbosity is a console filter.

## The safety we trade, and how we recover it

The one thing structs check that json doesn't is **key/type typos**.  Recover it:
1. Each section's keys live in ONE thin **writer helper** (`AddLadderLevel(json&, level, N, ecut, nG, role)`) —
   typed params, not stringly-typed at every call site — plus a **schema-check unit test** asserting the
   produced report's keys/shape.
2. **Units / precision** (needed regardless of the model): a small **format-hints** table (field name → unit +
   precision) the renderer consults; unknown fields get a `%g` default.  Units stay in the schema doc, never
   baked per value.

Both backends are already vendored: `submodules/tabulate` and `submodules/json` (nlohmann, already used in
`src/Pseudopotential/Imp/GTH_Potentials.C`).

## Vocabulary

- **Data providers** — the subsystems that PUT data in: the basis set, the `LASolver` (eigen / conditioning),
  the `GPW_Evaluator` (grids, caches), `SCFParams`.  A provider owns WHEN and WHAT it reports; it writes
  key-free into `CurrentRunReport`.
- **The sink** — `qchem.Reporting`: `GlobalReport` (all reports, keyed) + `CurrentRunReport` (the active run's
  report; the key-free write target).  Owns HOW data is keyed.
- **The renderers** — `RenderConsole` (terminal; generic table/tree) + `RenderJson` (`dump()`).  Own HOW data
  is laid out (inferred from structure).

## Ownership: a GLOBAL sink (like Cache and logging), NOT a threaded object

The dynamic-flow walkthrough settled this.  The actors:

| Actor | Class | Provides |
|---|---|---|
| 1 orchestrator  | `Calculation` (prod) / `RunGPW` (test) — builds the basis AND runs the SCF | owns lifetime; opens/keys/renders |
| 2 basis + grids | `GPWFactory` → `GPW_Evaluator` (built BEFORE the iterator) | `grids`, `basis.meta` |
| 3 eigen / cond  | `LASolver`, via `SCFIterator → WaveFunction::Factory → tCompositeWF → MakeIrrepWFs → SetBasisOverlap` | `basis.perIrrep`, `basis.removed` |
| 4 caches        | `GPW_Evaluator` stream/integral caches, during SCF | `cache` |
| 5 settings      | `SCFParams` (orchestrator holds it) | `scf` |

Actor 3 is the tell: the `LASolver` sits FIVE layers below the orchestrator, and none of the intervening
constructors (`WaveFunction::Factory`, `tCompositeWF`, `MakeIrrepWFs`) has any other reason to carry a
`report&`.  PUSHING a report parameter through them is signature pollution; PULLING needs getters wired up the
same five layers.  This is a cross-cutting concern, and the codebase already answers it TWICE the same way:
the integral **Cache is a singleton**, and the **`Report*` toggles** (`ReportOverlapConditioning`,
`ReportGridCharge`, `ReportBandGap`) are process-wide.  So `qchem.Reporting` is a **global sink**, initialised
once at app startup (singleton).

**Why `qchem.Reporting`, not `RunReport`:** `GlobalReport` holds more than SCF runs.  Example: the unit-test
suite does N runs, then dumps cross-run cache-reuse% stats — a report tied to NO single run.  So the module
name must not imply "run"; `GlobalReport` / `CurrentRunReport` are entities *within* `qchem.Reporting`.

**Providers write key-free; the orchestrator owns the key.**  A data provider writes at its CURRENT print site
(a one-line `cout` → `CurrentRunReport().EmitX(...)` swap) — NO new parameters, and it never sees a run key.
The orchestrator opens a run (`Begin(name)` → key `"<name>@<startTime>"`) and FILES `CurrentRunReport` into
`GlobalReport` under that key.  This is what dissolves the "each provider needs the key" problem.

**Nested runs (the HF/DHF SAD bootstrap, `project_numericcd_refactor`) don't clobber.**  `Begin`/`End` nest:
each run gets its own `CurrentRunReport` context and its own top-level `GlobalReport` key, and providers always
write to "the current" one.  The JSON KEEPS every run (parent + nested + sequential) for the GUI/HDF5; a depth
counter keeps the CONSOLE quiet for nested runs (render only at depth 1) — full record in JSON, tidy terminal.
Thread-safety is a non-issue: setup is single-threaded; OMP only parallelises the pair loops, which never
write the report.

**Incremental rendering — each section to the console AS IT COMPLETES, not batch-at-end.**  `EmitX` appends to
the current run's json AND (at console depth 1) immediately renders THAT section.  Two reasons: (1) BAIL-OUT
resilience — a crash mid-setup still shows every section that completed; (2) SLOW-section responsiveness —
basis/eigen appear at once instead of blocking on a slow grid build.  It costs nothing: the generic renderer
takes ANY json subtree, so "render as it lands" is just calling it on the section subtree.  The render call
needs only the **section name** — the current run is implicit (`CurrentRunReport`); the run key stays with the
orchestrator.

## Section schema (the JSON contract — freeze this before the GUI depends on it)

Top-level keyed by run (`"<name>@<startTime>"`) so nested + sequential runs coexist; within each run, ordered
sections; every field is a scalar / string / list thereof (JSON-native).  Shape each section for its inferred
layout (arrays-of-uniform-objects → tables; nested → trees):

```
{ "<name>@<startTime>" : {          // one entry per run (Begin/End open/close it)
  meta        { title, name, startTime, structureName, nElectrons, spinPolarized }
  basis       { name, engine, angular, nFunctions,
                perIrrep : [ { irrep, nFunctions, lambdaMin, lambdaMax, cond } ],
                removed  : [ { index, L, alpha, atom, position } ] }   // from PivotedCholeskyDrops (9b546bc1)
  grids       { densityEcut, cutoffFactor, raster,
                ladder : [ { level, N:[nx,ny,nz], ecut, nG, role } ],  // role: reference | coarse | rung
                localPP: { kappa, fieldSharpness } }
  cache       { tiers : [ { name, ramMB, pct } ],
                reuse : [ { cache, entries, lookups, reusePct, ramMB } ] }
  scf         { standard: { nMaxIter, minDrho, minDFD, minDE, minFD, mixer, kerkerG0,
                            accelerator, ortho, smearingkT },
                advanced: { pulayDepth, pulayStart, useMOM, momStartIter, momSmearPenalty,
                            guard:{holePersistence,maxReleases} } }
} }
```

Notes:
- `basis.removed` IS the "functions removed / fix-your-basis" report — named by exponent+angular+atom so the
  message is a two-minute basis-file edit, not "index 40".
- `scf.standard` vs `scf.advanced` mirrors the GPWPlan1 §1 Standard/Advanced grading (settled knobs vs
  pathological-case dials), so the display teaches the recipe.
- Keep field NAMES stable once the GUI consumes them — additive changes only (new fields, never renames).
- FUTURE sections: structure, symmetry group, irreps, and the Hamiltonian (also displayed to the console).

## Data providers → sections

| Section | Provided by | Replaces the ad-hoc print |
|---|---|---|
| `basis` | the basis + `qchem::PivotedCholeskyDrops` | `[overlap S]`, `[ortho] Auto…` |
| `grids` | `GPW_Evaluator` (ladder build) | `[GPW grid] ladder L=…` |
| `cache` | the evaluator stream/integral caches | `IntegralsCache RAM usage report` |
| `scf`   | `SCFParams` + mixer/accel `Tag()` | (nothing today) |

The writer helpers live WITH the `qchem.Reporting` leaf (not the providers): a provider passes its typed data
to `CurrentRunReport().EmitX(...)`; the provider owns WHEN/WHAT, the leaf owns HOW it's keyed.  All key names +
shapes are then greppable in one module next to the schema-check test.

This "run header" is distinct from the two displays that already exist / are planned:
1. **Run header** — this doc (basis, grids, cache, settings).  Rendered incrementally as setup proceeds.
2. **Per-iteration table** — GPWPlan1 §2, DONE (`DisplayColumns`/`DisplayColumnHeaders`): convergence per step.
3. **Final energy breakdown** — `EnergyBreakdown::Display` at run end.
Longer term all three are sections of one run *document* (the natural HDF5 / JSON record the GUI stores).

## Module / library placement

- `qchem.Reporting` — new leaf module.  Home: `src/Common/Reporting.C` (+ `Imp/`), beside the other
  cross-cutting leaves (`Types`, `Streamable`).  Interface = the `json` alias, `GlobalReport`/`CurrentRunReport`
  accessors, the per-section writer helpers, `RenderConsole`/`RenderJson`, and the format-hints table.
- tabulate + nlohmann are header-only; both already build in-tree — wire their includes into the `Reporting`
  target ONLY (keep them out of everything else).
- EVERY subsystem depends on `qchem.Reporting`, never the reverse — so no cycles (the compiler enforces this).

## Migration plan (incremental, one section per step, each independently shippable)

0. **Skeleton**: `qchem.Reporting` module — the `json` alias, `GlobalReport`/`CurrentRunReport` with the
   Begin/End key + depth machinery, the ONE generic `RenderConsole` (infer table/tree by depth+breadth, skip
   empty sections), `RenderJson`=`dump()`, and the format-hints table.  Unit test: hand-build a report json,
   render to console (golden-ish) and assert the schema (keys/shape) — the test IS the compile-time-check
   replacement.
1. **`scf`** (easiest — pure `SCFParams`): fill + render.  Proves the pipeline end-to-end; adds the
   currently-missing "what recipe am I running" header.
2. **`basis`**: fill from the basis + `PivotedCholeskyDrops`; retire `report_conditioning`'s line and the
   `[ortho]` warning into `basis.perIrrep` + `basis.removed`.  **Delivers the fix-your-basis report.**
3. **`grids`**: move the `[GPW grid]` ladder prints into `grids.ladder`.
4. **`cache`**: move the `IntegralsCache RAM usage report` into `cache`.
Each step deletes scattered `cout`s and adds one tidy block; behaviour is otherwise unchanged.

## JSON-for-GUI decisions to lock now

- **Stable names, additive evolution.**  Once the GUI reads a field, never rename/retype it; only add.
- **One document, versioned.**  Add a top-level `schemaVersion` so the GUI can guard.
- **Units explicit.**  Energies Hartree, lengths Bohr, exponents a.u. — state them in the schema doc; do NOT
  bake unit strings into every value (the renderer formats).
- **Render selection.**  The caller (CLI vs binding) picks `RenderConsole` vs `RenderJson`; the model is identical.

## Open questions / styling

- Tabulate styling: terminal-portable default, colour opt-in; professional, not busy.
  - When screen real estate is tight (e.g. the per-iteration table later) drop the row/column dividers.
  - Use UTF liberally for Greek symbols etc.
  - For atoms, colour distinguished the irreps (s,p,d,f) and looked great — do similar for polarised Gaussians.
- FUTURE: does the run header ALSO belong in a per-system `DisplayColumnHeaders`-style virtual (Molecular vs
  Solid add/drop sections), or is it a flat top-level assembly?  Lean: flat assembly (sections are universal;
  only their CONTENT varies), empty sections skipped (no grids on the molecular path).
- FUTURE: HDF5 — is the JSON report a sidecar, or serialised INTO the HDF5 run group?  Defer; `RenderJson`
  makes either trivial.

## Pointers
- Detector already committed: `qchem::PivotedCholeskyDrops` (`9b546bc1`) — feeds `basis.removed`.
- Per-iteration display (the sibling report): `doc/GPWPlan1.md` §2 (DONE).
- GUI/JSON consumer: `project_viz_gui_plan` (nanobind + PySide6 + pyqtgraph + HDF5).
- Standard/Advanced knob grading: `doc/GPWPlan1.md` §1.
