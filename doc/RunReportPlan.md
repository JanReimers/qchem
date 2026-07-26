# RunReportPlan — a professional, JSON-ready run report

Status: PLAN (2026-07-25).  Consolidate the scattered setup diagnostics into one organized "run header",
rendered as tidy terminal tables (tabulate) today and JSON for the GUI tomorrow — from a single data model.

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

## The core principle: the report IS json; ONE generic renderer

Separate the report's **data** from its **rendering** — and take the data model all the way to
`nlohmann::json` itself (user, 2026-07-26).  No section structs, no `to_json`/`from_json`, no boilerplate:

```
RunReport  ==  an nlohmann::json object          (a LEAF module — depends only on json + Types)
     |         written via thin per-section WRITER HELPERS (typed args; the key strings live in ONE place)
     |
     +--> RenderText(const json&, ostream&)   -> tabulate  (terminal, today; ONE generic walker)
     +--> RenderJson(const json&)             -> report.dump()  (GUI, tomorrow; trivial)
```

**Why json-as-model, not structs:** a report is DISPLAY, not physics.  A wrong key is a cosmetic bug, not a
wrong number — so the project's "prefer compile-time" bias (which governs the wavefunction types) is the wrong
tool here.  The line we hold: **typed physics, schemaless display.**  Payoff: JSON output is literally
`report.dump()`; adding a field is writing a key; the GUI contract is the json directly.

**The text renderer is NOT an extra cost of json — you write table-building code either way** (user,
2026-07-26).  Structs wouldn't save it: C++ has no reflection, so with hard-coded types you hand-write a
renderer PER section and add one for every new section.  json is runtime-introspectable, so a generic walker
can render the sections that obey the shape conventions below — the same introspectability that "loses"
compile-time key checking is what buys the generic path.

**But temper the optimism (user, 2026-07-26): a generic json→pretty-table walker is NOT a universal magic
printer.**  It handles the CLEAN shapes well; it does NOT auto-produce professional output for arbitrary json
(non-uniform arrays, nested tables, wide-for-the-terminal tables, values wanting special formatting).  Two
things make it work in practice, and neither is "the renderer is clever":
- **Design the json layout WITH TABLES IN MIND** — deliberately make each section an array-of-uniform-objects
  or a flat scalar map, so it FALLS INTO a shape the walker prints well.  The layout does the work.
- **Per-section render OVERRIDE as an escape hatch** — a section that genuinely doesn't fit registers its own
  small renderer.  So: generic by DEFAULT, override where needed — still far less code than hand-rendering
  everything, without pretending the generic path covers every case.  "We shall see" as the sections land.

The two shapes the default walker keys off:
- **array of uniform objects → a tabulate table** (object keys = column headers, elements = rows):
  `grids.ladder`, `basis.perIrrep`, `basis.removed`, `cache.reuse`.
- **object of scalars → a two-column key/value table**: `scf.standard`, `meta`.
- nesting → a subsection header, then recurse.

**The ONE thing structs actually check that json doesn't — key/type typos** — plus the formatting every model
needs regardless:
1. *Key typos* (the genuine json-specific gap): each section's keys live in ONE thin writer helper
   (`AddLadderLevel(json&, level, N, ecut, nG, role)`), typed params, not scattered stringly at every call
   site; plus a schema-check unit test asserting the produced report's keys/shape.  ~All the safety, no struct.
2. *Units / precision* (needed with structs too — a per-section renderer would hard-code it): a small
   **format-hints** table (field name → unit + precision) the generic renderer consults; unknown fields get a
   `%g` default.  Units stay in the schema doc, never baked per value.

Both backends are already vendored: `submodules/tabulate` and `submodules/json` (nlohmann, already used in
`src/Pseudopotential/Imp/GTH_Potentials.C`).

## Section schema (the JSON contract — freeze this before the GUI depends on it)

Ordered sections; every field is a scalar / string / list of scalars (JSON-native):

```
RunReport
  meta        { title, timestamp?, structureName, nElectrons, spinPolarized }
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
                            guard:{holePersistence,maxReleases}, momSmearPenalty, … } }
```

Notes:
- `basis.removed` IS the "functions removed / fix-your-basis" report — it lands here, in its proper home,
  named by exponent+angular+atom so the message is a two-minute basis-file edit, not "index 40".
- `scf.standard` vs `scf.advanced` mirrors the GPWPlan1 §1 Standard/Advanced grading (settled knobs vs
  pathological-case dials), so the display teaches the recipe.
- Keep field NAMES stable once the GUI consumes them — additive changes only (new fields, never renames).

## Architecture: a leaf, written by everyone via helpers

`qchem.RunReport` is a **leaf** — the `json` typedef alias, the per-section writer helpers, the two renderers,
and the format-hints table; depends only on `qchem.Types` + tabulate/json.  EVERY subsystem depends on it,
never the reverse — so no cycles (the compiler enforces this; see CLAUDE.md).  Each subsystem calls ONLY its
own section's writer helper (which owns that section's key names):

| Section | Writer called by | Replaces the ad-hoc print |
|---|---|---|
| `basis`  | the basis + `qchem::PivotedCholeskyDrops` | `[overlap S]`, `[ortho] Auto…` |
| `grids`  | `GPW_Evaluator` (ladder build) | `[GPW grid] ladder L=…` |
| `cache`  | the evaluator stream/integral caches | `IntegralsCache RAM usage report` |
| `scf`    | `SCFParams` + mixer/accel `Tag()` | (nothing today) |

Keeping the writers WITH the report leaf (not the subsystems) means the key strings + shapes live in one
module — the schema is greppable in one file, and the schema-check test sits right next to them.

The `Calculation` / `tSCFIterator` assembles the report at run start and renders it ONCE — the "run header".
It is distinct from the two displays that already exist or are planned:
1. **Run header** — this doc (setup: basis, grids, cache, settings).  Rendered once, at start.
2. **Per-iteration table** — GPWPlan1 §2, DONE (`DisplayColumns`/`DisplayColumnHeaders`): convergence per step.
3. **Final energy breakdown** — the `EnergyBreakdown::Display` at run end.
Longer term all three are sections of one run *document* (also the natural HDF5 / JSON record the GUI stores);
for now, ship the header.

## Module / library placement

- `qchem.RunReport` — new leaf module.  Suggested home: `src/Common/RunReport.C` (+ `Imp/`), beside the other
  cross-cutting leaves (`Types`, `Streamable`).  Interface = the `json` alias, the per-section writer helpers,
  `RenderText`/`RenderJson`, and the format-hints table.
- tabulate is header-only; nlohmann is header-only.  Both already build in-tree — wire them into the
  `RunReport` target's includes only (keep them out of everything else).
- The writer helpers live WITH the report leaf (not the subsystems): a subsystem passes its typed data to
  `report::AddGridLevel(...)` etc., so all key names + shapes are greppable in one file next to the schema
  test.  The subsystem still owns WHEN/WHAT to report; the leaf owns HOW it's keyed.

## Migration plan (incremental, one section per step, each independently shippable)

0. **Skeleton**: `qchem.RunReport` module — the `json` alias, the ONE generic `RenderText` walker (array-of-
   objects → table, object-of-scalars → key/value, recurse; skip empty sections), `RenderJson`=`dump()`, and
   the format-hints table.  Unit test: hand-build a report json, render to text (golden-ish) and assert the
   schema (keys/shape) — the test IS the compile-time-check replacement.
1. **`scf`** (easiest — pure `SCFParams`): fill + render.  Proves the pipeline end-to-end; adds the
   currently-missing "what recipe am I running" header.
2. **`basis`**: fill from the basis + `PivotedCholeskyDrops`; retire `report_conditioning`'s ad-hoc line and
   the `[ortho]` warning into `basis.perIrrep` + `basis.removed`.  **Delivers the fix-your-basis report.**
3. **`grids`**: move the `[GPW grid]` ladder prints into `grids.ladder`.
4. **`cache`**: move the `IntegralsCache RAM usage report` into `cache`.
Each step deletes scattered `cout`s and adds one tidy tabulate block; behaviour is otherwise unchanged.

## JSON-for-GUI decisions to lock now

- **Stable names, additive evolution.**  Once the GUI reads a field, never rename/retype it; only add.
- **One document, versioned.**  Add a top-level `schemaVersion` so the GUI can guard.  (Cheap now, priceless
  later.)
- **Units explicit.**  Energies Hartree, lengths Bohr, exponents a.u. — state them in the schema doc; do NOT
  bake unit strings into every value (the GUI formats).
- **Render selection.**  A single `enum Format { Text, Json }` (or two `Render` overloads) — the caller (CLI
  vs binding) picks; the model is identical.

## Open questions

- tabulate styling: how much (row separators, colour, unicode box-drawing)?  Terminal-portable default,
  colour opt-in.  Keep it professional, not busy.
- Does the run header belong ALSO in the per-system `DisplayColumnHeaders` virtual (so Molecular vs Solid can
  add/drop sections), or is it a flat top-level assembly?  Lean: flat assembly (the sections are universal;
  only their CONTENT varies), with empty sections skipped (no grids on the molecular path).
- HDF5: the GUI plan stores runs in HDF5 — is the JSON report a sidecar, or serialized INTO the HDF5 run
  group?  Defer; the JSON renderer makes either trivial.

## Pointers
- Detector already committed: `qchem::PivotedCholeskyDrops` (`9b546bc1`) — feeds `basis.removed`.
- Per-iteration display (the sibling report): `doc/GPWPlan1.md` §2 (DONE).
- GUI/JSON consumer: `project_viz_gui_plan` (nanobind + PySide6 + pyqtgraph + HDF5).
- Standard/Advanced knob grading: `doc/GPWPlan1.md` §1.
