# RunReportPlan — a professional, JSON-ready reporting layer

Status: SKELETON + ALL FOUR SECTIONS (`scf`, `basis`, `grids`, `cache`) DONE (2026-07-26) — the migration is
COMPLETE.  Consolidate the scattered setup diagnostics into one organized, machine-readable report — rendered as
tidy terminal output today and JSON for the GUI tomorrow, from a single data model.

Post-migration polish (2026-07-26): unused cache tiers (0 RAM) and untouched reuse caches are dropped; the legacy
dtor `percent(0,0)` `-nan%` guarded to `0%`; `FormatHint.fixed` gives clean fixed-point percentages/MB.  The GPW
`[stream cache]` collocation-coverage line is folded into `grids.stream` via a new `report::EmitAt(section, key,
value)` (write a late sub-block into an already-emitted section + render just that block; it builds after the
grids table renders).  Remaining/future work: the SOLID refactor (deferred, see "SOLID target"),
`basis.removed` exponent/atom naming, disk/rolling-log sinks, and a NON-reporting finding the report surfaced —
the SCF setup builds grids BEFORE the overlap-conditioning eigen-analysis (a potential show-stopper), so a
singular basis is detected only after the grid ladder is built; a fail-fast reorder is worth considering.

**Step 0 landed**: `qchem.Reporting` leaf module in `qcCommon` (`src/Common/Reporting.C` + `Imp/Reporting.C`) —
`json = nlohmann::ordered_json` (ordered so sections keep emit order), global sink (`GlobalReport` keyed +
key-free `CurrentRunReport`), `Begin`/`End` nesting with a depth counter, `SetConsole`/`EmitSection` (incremental
depth-1-only console render), the ONE generic `RenderConsole` (layout inferred: key/value table | multi-column
table | indented tree), `RenderJson=dump(2)`, and the `HintFor` format-hints table.  Schema-check + render tests
in `src/Common/tests/Reporting.C` (9 green, part of `UTCommon`), plus a `DISABLED_VisualDump` for eyeballing.

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

### Renderer vs SINK — where disk / rolling logs fit (2026-07-26)

The "third renderer type for (rolling?) log files on disk" is best NOT modelled as a third renderer.  Two
orthogonal axes:
- **Renderer** = *json → representation*: `RenderConsole` (text, table/tree) and `RenderJson` (`dump()`).
- **Sink** = *where the representation goes*: stdout, a file, a rolling file, an HDF5 group.

A disk log is **(a renderer) × (a file sink)**, so `RenderConsole(json, ostream&)` already takes the sink as an
`ostream&` — point it at an `ofstream` and you have a text sidecar; `RenderJson` gives the `report.json` sidecar
for free.  Two disk stories, mapping onto different entities:
1. **Per-run sidecar** (`report.json` / `report.txt`) — from `CurrentRunReport` at `End()`.  Essentially free.
2. **Rolling cross-run log** — append-only, time-ordered, spanning runs, rotated by size/date.  Lives on
   `GlobalReport`, not `CurrentRunReport`.  The *rotation* logic (open → size-check → rotate → retain) is a small
   `RollingFileSink` that consumes whatever a renderer emits — it must NOT leak into `RenderConsole`.

Incremental section-by-section rendering composes directly with a file sink (tee the same stream to disk; the
crash-resilience argument gets stronger on disk).  So the seam is already right; **rotation policy is DEFERRED**
until a concrete need — it drops in later as a `GlobalReport` sink, no renderer change.

### SOLID target (DEFERRED — refactor after migration steps 1-4 land, when real usage is known)

The skeleton is deliberately a flat namespace of free functions to get sections migrated and *see what lands*
before committing to abstractions.  The target once the provider/client usage is concrete:
- **ISP — two narrow facets over the one global sink.**  A *provider* facet (`CurrentRunReport()` + the typed
  `EmitX(...)` writer helpers) is ALL a data provider sees — never Begin/End, GlobalReport, or renderers.  An
  *orchestrator/client* facet (Begin/End lifecycle + render selection) is what `Calculation`/`RunGPW` see.
  Same backing singleton, two segregated interfaces.
- **DIP/LSP — an abstract `Renderer`.**  `class Renderer { virtual void Render(const json&, std::ostream&) = 0; };`
  with concrete `ConsoleRenderer` / `JsonRenderer` (+ future `Html`/`Xml`/`Binary`/`Log`), selected by a factory
  `std::unique_ptr<Renderer> MakeRenderer(RendererType{Console,Json,Html,Log,Xml,Binary,...})`.  Client code
  depends on the abstraction, not the concretes.  `RenderConsole`/`RenderJson` become the method bodies.
- **Renderer vs sink stays orthogonal here too.**  Console/Json/Html/Xml/Binary/`Log` are serialization FORMATS
  (renderers: json→bytes); the `ostream&` is the SINK (stdout / file / rolling).  `Log` is a *line-oriented
  renderer* (one timestamped line per section); rotation lives in the sink (`RollingFileSink`), not in `Log`.
- **Why defer:** with only Console+Json real today, a hierarchy + factory would be speculative (YAGNI/`HTML`
  et al. have no concrete requirements yet); the call sites are few and internal, so the lift is cheap later.

`Begin(name)` auto-stamps `NowIso()` (production ergonomics); `Begin(name, startTime)` is the deterministic
test/repro seam.

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
1. **A per-section schema-check test.**  The data provider builds the section json DIRECTLY (the report IS json —
   no mediating C++ display struct; a struct is just a parallel field-list that duplicates the provider's own
   state and drags domain vocabulary into the `qchem.Reporting` leaf).  The provider owns the key strings, right
   where the section's semantics are understood; a **schema-check test** on the emitted json is the typo-catcher
   that stands in for the compile-time check.  (DECISION 2026-07-26: an earlier draft routed each section through
   a typed `EmitX(struct)` writer in the leaf — dropped.  The provider facet is just `{CurrentRunReport(),
   EmitSection(name, json)}`, fully generic; this is the same compile-time→runtime compromise the codebase
   already makes with cross-interface `dynamic_cast`.)
2. **Units / precision** (needed regardless of the model): a small **format-hints** table (field name → unit +
   precision) the renderer consults; unknown fields get a `%g` default.  This is the ONE place domain field
   NAMES appear on the render side, and it is display config, not data — units stay here + in the schema doc,
   never baked per value.

Both backends are already vendored: `submodules/tabulate` and `submodules/json` (nlohmann, already used in
`src/Pseudopotential/Imp/GTH_Potentials.C`).

## The section cursor — how a five-layers-down provider writes with NO threading (added 2026-07-26, step 2)

`Begin`/`End` lets a deep provider write to "the current RUN" implicitly (it never sees a run key).  Step 2's
`basis` section needs the same trick one level finer: the conditioning numbers are computed in the `LASolver`,
five layers below the orchestrator (`SCFIterator → WaveFunction::Factory → tCompositeWF → MakeIrrepWFs →
SetBasisOverlap`), and the `LASolver` does not know which irrep it is solving.  Rather than thread a `report&`
or an irrep label down through every intervening constructor — signature pollution the global sink exists to
avoid, and a problem that recurs at EVERY layer — the **sink holds a navigation cursor**:

- The sink stores WHERE the next write lands (a **path** of object-key / array-index segments, NOT a raw
  `json*` — `ordered_json` is vector-backed, so an insert can reallocate and dangle a pointer; we re-resolve
  from the run root each write, depth ≤ 3).
- An **ancestor** that knows the context opens it with an RAII handle: `report::Section basis("basis")` (a
  top-level section, renders once on close), `report::Row row("perIrrep")` (append an array element, descend).
- **Deep code writes context-free**: `report::Set("cond", x)` lands in the current run AND the current irrep
  row, automatically.  `MakeIrrepWFs` opens the `Row`; `LASolver::report_conditioning` just `Set`s into it.
- All ops are **inert when no run is open** (`Depth()==0`), so a kernel (`LASolver`) used outside a report —
  every LASolver unit test, `BandStructure`, `A_Pool` — is completely unaffected.
- RAII (not free Push/Pop) so an exception mid-assembly cannot leave the cursor dangling; the handles live ONLY
  in the context-openers, so the deep provider still gets zero new parameters.  Nested runs (SAD bootstrap)
  save/restore the outer cursor at `Begin`/`End`.

This is the mechanism the "global sink" was always for; step 2 is where it earned its keep.  Cursor API:
`Set(key, value)` + the `Section` / `Row` RAII handles (`src/Common/Reporting.C`).  Use the cursor for a section
ASSEMBLED across layers (`basis`); use `EmitSection(name, json)` for one a single provider builds whole (`scf`).

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

The provider builds its section json directly and calls `CurrentRunReport()`/`EmitSection(name, json)`: the
provider owns WHEN/WHAT **and the key strings** (co-located with the section's semantics), the leaf owns HOW it's
filed + rendered.  Keys are greppable at the provider; the per-section schema-check test freezes the shape.
(Superseded the earlier "typed `EmitX` writer helper in the leaf" idea — see "The safety we trade".)

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

0. **Skeleton** — ✅ DONE (2026-07-26).  `qchem.Reporting` module — the `json` alias (`ordered_json`),
   `GlobalReport`/`CurrentRunReport` with the Begin/End key + depth machinery, `SetConsole`/`EmitSection`
   (incremental depth-1 render), the ONE generic `RenderConsole` (infer table/tree by depth+breadth, skip empty
   sections), `RenderJson`=`dump()`, and the format-hints table.  Unit test hand-builds a report json, renders to
   console and asserts the schema (keys/shape) — the test IS the compile-time-check replacement.  NOTE:
   `EmitSection` is the generic key-free emit; the typed `EmitX` writer helpers land WITH each section below.
1. **`scf`** — ✅ DONE (2026-07-26).  The PROVIDER-FACET pattern in the flesh, and it is MINIMAL: the provider
   (`Calculation::Converge`, via the `EmitScfSection` static) builds the `scf` section json DIRECTLY from
   `SCFParams` + the accelerator tag (mixer DERIVED: Pulay>Kerker>linear) and calls `report::EmitSection("scf",
   scf)` — no display struct, no typed leaf writer (that draft was dropped; see "The safety we trade").
   `Converge` brackets the run (`Begin(structure.ID())` + an RAII `End()` guard so a throwing `Iterate` can't
   leak the run-context stack).  The facade only RECORDS (into `GlobalReport`); `scfrun` opts into console
   rendering via `report::SetConsole(cout)` — the first dogfood.  Schema-check test = `M_Calculation.
   ScfReportSchema` (with the provider, asserts the emitted keys/shape); full `ctest -j16` green (599/599).
   NOTES: (a) `ortho`/`accelerator` are additive; `ortho` fills when the `basis`/LASolver section lands.
   (b) `minDE=1e30` sentinel renders as `1e+30` ("off") — a future display nicety could show "off".  (c) scfrun
   shows the header twice (facade ctor auto-converges, then scfrun re-converges) — each `Converge` is correctly
   its own run, not a bug.
2. **`basis`** — ✅ DONE (2026-07-26).  The CURSOR path (see "The section cursor" above).  `MakeIrrepWFs` opens
   a `report::Row("perIrrep")` per irrep block and stamps `irrep`/`nFunctions`; `LASolver::report_conditioning`
   `Set`s `lambdaMin`/`lambdaMax`/`cond` into that row from five layers down (its `[overlap S]` console line
   kept behind the existing default-OFF toggle, now feeding the report whenever a run is open).  `basis.removed`
   comes from `PivotedCholeskyDrops` per block (empty on a healthy basis).  The orchestrator opens the
   `Section("basis")` and stamps `name`/`engine`/`angular`/`nFunctions`.  Schema test `M_Calculation.
   BasisReportSchema` (symmetry ON → multi-row `perIrrep` table); `ctest -j16` green (604/604, timing flat).
   LIMITATION: `basis.removed` names by `{irrep, index}` only — the schema's `L`/`alpha`/`atom`/`position`
   naming ("the fix-your-basis report") needs a per-function metadata accessor on the block that does not exist
   yet; that enrichment is a §4a-actuator follow-up.  The `[ortho]` drop warnings still `cerr` (not yet retired;
   `removed` now carries the same information structurally).
3. **`grids`** — ✅ DONE (2026-07-26).  `GPW_Evaluator::EmitGridsReport` builds the whole section
   (`densityEcut`/`cutoffFactor`/`raster` + the `ladder` rows {`level`,`N`,`ecut`,`nG`,`role`=reference|coarse|
   rung} + `localPP.kappa`) and `EmitSection`s it — single-provider, so NO cursor.  The basis-ctor call site
   (`Lattice_3D/Imp/BasisSet.C`) routes to `EmitGridsReport` when a run is open, else the legacy
   `ReportGrids(cout)` (unbracketed basis-only tests unchanged).  **This required bracketing the GPW run** —
   GPW does not use the molecular `Calculation` facade, so `RunGPW` (the GPW orchestrator, per the ownership
   table) now `Begin`/`End`s the run and opens `Section("basis")` around the `SolidSCFIterator` construction.
   BONUS: that Section means GPW ALSO gets `basis.perIrrep` — per-Bloch-block conditioning through the SAME
   cursor path as step 2 (the `LASolver` `dcmplx` write).  Schema test `GPW_SCF.GridsReportSchema`; `ctest -j16`
   green (605/605).
4. **`cache`** — ✅ DONE (2026-07-26).  PER-RUN SNAPSHOT (user decision): the `IntegralsCache` is a process-wide
   singleton whose `~dtor` `ReportRAMUsage(cout)` is CUMULATIVE across all runs, not per-run.  So a new
   `IntegralsCache<T>::EmitReport()` (virtual, built in `IntegralsCache_RAM`) snapshots the current stats into
   the open run's `cache` section (`tiers` [Jac/Kab/Cach4 ramMB+pct] + `reuse` [per Cache2/3/4:
   entries/lookups/reusePct/ramMB]); `Calculation::Converge` calls it AFTER `Iterate`, inside the bracket.  The
   cache is NEVER cleared between runs (user pin), so each snapshot is cumulative-to-that-point — exact for a
   one-run process, the running total otherwise; the legacy dtor `cout` stays as the process-level summary.
   Added `FormatHint.fixed` (decimal places vs sig-figs) so percentages/MB render `50.0`/`99.9`/`28.078`, not
   `5e+01`.  Schema test `M_Calculation.CacheReportSchema`; `ctest -j16` green (606/606).  NOTE: this is the
   MOLECULAR cache (`theCache<double>`); the GPW `[stream cache]` is a separate cache, a future section.
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
