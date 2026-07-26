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

## The core principle: DATA MODEL, not print statements

The single decision that makes this clean — and makes JSON *free* rather than a parallel code path — is to
separate the report's **data** from its **rendering**:

```
RunReport            (pure data; a LEAF module — depends on nothing but Types)
   Section structs   (numbers + strings + vectors; NEVER pre-formatted lines)
        |
        +--> Render(report, ostream)  -> tabulate tables + section rules   (terminal, today)
        +--> Render(report, json&)    -> nlohmann::json                     (GUI, tomorrow)
```

Because sections hold values, not strings, the terminal and the GUI render the *identical* model — they can
never drift.  This is why we design for JSON now even though only the text renderer ships first: the schema is
the contract, and getting it right up front is the whole point.

Both renderers' backends are already vendored: `submodules/tabulate` and `submodules/json` (nlohmann, already
used in `src/Pseudopotential/Imp/GTH_Potentials.C`).

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

## Architecture: a leaf, filled by everyone

`RunReport` is a **leaf** (data + the two renderers; depends only on `qchem.Types` + tabulate/json).  EVERY
subsystem depends on it, never the reverse — so no cycles (the compiler enforces this; see CLAUDE.md).  Each
subsystem fills ONLY its own section:

| Section | Filled by | Replaces the ad-hoc print |
|---|---|---|
| `basis`  | the basis + `qchem::PivotedCholeskyDrops` | `[overlap S]`, `[ortho] Auto…` |
| `grids`  | `GPW_Evaluator` (ladder build) | `[GPW grid] ladder L=…` |
| `cache`  | the evaluator stream/integral caches | `IntegralsCache RAM usage report` |
| `scf`    | `SCFParams` + mixer/accel `Tag()` | (nothing today) |

The `Calculation` / `tSCFIterator` assembles the report at run start and renders it ONCE — the "run header".
It is distinct from the two displays that already exist or are planned:
1. **Run header** — this doc (setup: basis, grids, cache, settings).  Rendered once, at start.
2. **Per-iteration table** — GPWPlan1 §2, DONE (`DisplayColumns`/`DisplayColumnHeaders`): convergence per step.
3. **Final energy breakdown** — the `EnergyBreakdown::Display` at run end.
Longer term all three are sections of one run *document* (also the natural HDF5 / JSON record the GUI stores);
for now, ship the header.

## Module / library placement

- `qchem.RunReport` — new leaf module.  Suggested home: `src/Common/RunReport.C` (+ `Imp/`), beside the other
  cross-cutting leaves (`Types`, `Streamable`).  Interface = the section structs + `Render` overloads.
- tabulate is header-only; nlohmann is header-only.  Both already build in-tree — wire them into the
  `RunReport` target's includes only (keep them out of everything else).
- The section-FILL helpers live with their subsystem (e.g. `GridSection FillGrids(const GPW_Evaluator&)`),
  so `RunReport` stays basis/grid-agnostic and each subsystem owns its own schema mapping.

## Migration plan (incremental, one section per step, each independently shippable)

0. **Skeleton**: `qchem.RunReport` module — the section structs + a text renderer (tabulate) that prints the
   sections it has, skipping empty ones + a JSON renderer stub (`to_json` per section).  Unit test: build a
   populated report, render to text (golden-ish) and to JSON (schema check).
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
