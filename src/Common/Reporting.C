// File: Common/Reporting.C  A JSON-ready reporting layer -- the report IS json.
//
//   Data providers write key-free into CurrentRunReport(); renderers turn the
//   json into console text (a single GENERIC renderer -- table vs tree, layout
//   INFERRED from the json's structure, NO per-section renderers) or a json
//   string for the GUI/HDF5.  See doc/RunReportPlan.md for the full design.
//
//   Typed physics, schemaless display: a report is display, not physics, so we
//   trade compile-time key-safety for an extensible schema and recover it with
//   thin typed writer helpers + a schema-check unit test.
module;
#include <string>
#include <cstddef>
#include <iosfwd>
#include <nlohmann/json.hpp>

export module qchem.Reporting;

export namespace qchem::report {

//! The report data model IS an nlohmann json.  ordered_json keeps sections in the
//! order providers emit them (a plain json object would sort keys alphabetically).
using json = nlohmann::ordered_json;

//! Console verbosity.  A CONSOLE-ONLY concern: providers ALWAYS store everything,
//! so the json/GUI/HDF5 record is always complete; the level only decides how much
//! of that complete report RenderConsole SHOWS.  (Orthogonal to layout, which comes
//! from structure.)
enum class Detail { Terse, Normal, Verbose };

//============================================================================
// The sink -- a GLOBAL singleton (like the integral Cache and the Report* toggles).
//
// Providers write key-free into CurrentRunReport() at their current print site --
// no new parameters, they never see a run key.  The orchestrator brackets a run
// with Begin/End: Begin pushes a fresh current report, End files it into
// GlobalReport() under "<name>@<startTime>".  Runs NEST (the HF/DHF SAD bootstrap):
// each Begin/End pair is its own context and its own top-level GlobalReport key, so
// nested + sequential runs coexist in the json; a depth counter keeps the CONSOLE
// quiet below depth 1.  Setup is single-threaded (OMP only parallelises pair loops,
// which never write the report), so no locking is needed.
//============================================================================

//! The active run's document -- the key-free write target for data providers.
//! Undefined behaviour to call outside a Begin/End bracket (asserts in debug).
json& CurrentRunReport();

//! Every run ever begun, keyed "<name>@<startTime>" (the GUI/HDF5 record).
const json& GlobalReport();

//! Drop all accumulated runs from GlobalReport (a test/utility hook, like Cache::Clear;
//! also the seam a long-lived process uses to bound retention).  No run may be open.
void ClearGlobal();

//! Open a run.  The common form auto-stamps startTime with NowIso() -- production
//! callers should NOT have to remember a timestamp.
void Begin(const std::string& name);

//! Explicit-timestamp overload: the deterministic seam for tests / reproducible
//! runs (a fixed string keeps a golden-ish render stable; the wall clock does not).
void Begin(const std::string& name, const std::string& startTime);

//! Close the current run: file it into GlobalReport under its key, pop the context.
void End();

//! Current nesting depth (0 == no run open, 1 == top-level run).  Console renders
//! only at depth 1; the json keeps every depth.
int Depth();

//! True if a section named `name` is currently open on the cursor (an ancestor opened
//! a Section/Row path through it).  Lets a deep provider emit ONLY when the orchestrator
//! has established the expected context -- e.g. conditioning writes into basis.perIrrep
//! only when a "basis" section is open, so the same code stays silent when it is not.
bool InSection(const std::string& name);

//! An ISO-8601-ish "YYYY-MM-DDThh:mm:ss" stamp for Begin's startTime (production use).
std::string NowIso();

//============================================================================
// Console configuration -- set ONCE at app startup, like a logging config.  Lets
// providers Emit key-free with no ostream parameter (the console is implicit).
//============================================================================

//! Point incremental console rendering at a stream + set the detail level.
void SetConsole(std::ostream& os, Detail = Detail::Normal);

//! A live PROGRESS heartbeat: flush a one-line "step starting" message to the console IMMEDIATELY (not stored
//! in the report json -- it is a pulse, not a record).  The structured sections render only when COMPLETE, so a
//! slow step (grid build, basis construction) is silent meanwhile; Log lets a long run show what it is doing
//! RIGHT NOW so the user can judge whether to let it run or kill + retune.  No-op when no console is attached.
//! (First increment of the deferred "logger" sink; a rolling-file logger is the natural next home for these.)
void Log(const std::string& message);

//! Detach the console (Emit still records into the json, just renders nothing).
void ClearConsole();

//============================================================================
// Timing ledger -- a process-global {label -> accumulated wall seconds} map, so ANY layer can charge
// its cost to a named bucket with no plumbing (the same no-refs-threaded philosophy as the cursor).
// The point is DISCOVERY: instrument every suspected setup/iteration site, emit once at run end sorted
// by cost, and whatever dominates simply shows up (the "where did my 300 s go" report).
//============================================================================

//! Accumulate \a seconds of wall time under \a key (creates the bucket on first use), and COUNT the call.
//! The count is what makes two buckets comparable: a total is only meaningful against another total when
//! both fire at the same cadence, and a per-iteration bucket beside a once-per-run one reads as a win that
//! is not there.  \c EmitTimings prints the per-call price for any bucket that fired more than once.
//! \note Entries charging no time (a scope whose body delegated entirely to a nested \c Timed) are NOT
//! counted -- the time is exclusive, so counting them would divide real work by inflated entries.
void AddTime(const std::string& key, double seconds);

//! RAII stopwatch: charges the scope's EXCLUSIVE wall time -- elapsed MINUS any enclosed Timed scopes --
//! to \a key on destruction (a thread-less stack tracks nesting).  So the emitted buckets are DISJOINT
//! and sum to the instrumented wall time: an enclosing phase ("seed + ortho") whose cost is really its
//! lazy children ("local-PP", "streams") shows only its own residue, and the sorted-by-cost table needs
//! no indentation to read truthfully.
class Timed
{
public:
    explicit Timed(std::string key);
    ~Timed();
    Timed(const Timed&)            = delete;
    Timed& operator=(const Timed&) = delete;
private:
    std::string itsKey;
    long long   itsT0;   //!< steady_clock ns at construction (opaque here; <chrono> stays in the impl)
};

//! Emit the ledger as top-level section \a name (labels sorted by cost, seconds) and CLEAR it.
//! A bucket entered more than once carries \c [xN, s/call] in its label -- the per-call price, which is
//! the number an optimisation decision needs (see \c AddTime).  No-op when the ledger is empty.
//! Call once at run end (after the SCF), when the lazy builds are done.
void EmitTimings(const std::string& name="timing");

//============================================================================
// Emit -- append a COMPLETED section to the current run AND (at depth 1, if a
// console is set) render THAT section immediately.  Incremental by design:
// (1) bail-out resilience -- a crash mid-setup still shows completed sections;
// (2) slow-section responsiveness -- fast sections appear without blocking on a
// slow grid build.
//============================================================================
//! Emit a completed section: a DATA PROVIDER builds the section json itself (the
//! report IS json -- no mediating C++ struct) and hands it here key-free.  Providers
//! own WHEN + WHAT (the keys); this leaf owns HOW it is filed + rendered.  The
//! provider facet is exactly {CurrentRunReport(), EmitSection()} -- fully generic,
//! no per-section methods, so `qchem.Reporting` carries no domain vocabulary.  Typo
//! safety is recovered by a per-section schema-check test on the emitted json.
void EmitSection(const std::string& name, json section);

//----------------------------------------------------------------------------
// Section cursor -- the SINK holds WHERE the next write lands, so a deep provider
// writes CONTEXT-FREE: no run key, no section path, no `report&` threaded through
// the physics code.  An ANCESTOR establishes context with an RAII Section/Row; code
// it calls into (however many layers down) just calls Set(key,value) and the field
// lands in the current run AND the current section/row automatically.  This is the
// Begin/End run-nesting idea taken one level finer -- and it is why a five-layers-
// down provider (LASolver) needs no parameters.  ALL cursor ops are inert when no
// run is open (Depth()==0), so a kernel used outside a report is unaffected.
//----------------------------------------------------------------------------

//! Write a field at the current cursor node (a no-op when no run is open).  `value`
//! is any json-convertible primitive (number/bool/string); stringify domain types at
//! the call site (the leaf carries no domain types).
void Set(const std::string& key, json value);

//! Write a sub-block `value` under `section[key]` and render JUST that sub-block to the console (heading
//! "section ▸ key").  For a late addendum to an ALREADY-emitted section -- e.g. a grids.stream that
//! builds after the grids table rendered.  A no-op when no run is open.  Cursor-independent (absolute path).
//! `minLevel` gates ONLY the console render (the json ALWAYS records the block): a Verbose-only sub-block
//! (e.g. basis.usage) stays in the record at every level but prints only when SetConsole's detail >= minLevel.
//! IDEMPOTENT: an identical value already at the path is a complete no-op (no rewrite, no re-render) -- so a
//! provider may announce unconditionally at its own trigger and the dedup lives here, RUN-scoped, instead of
//! in process-global latches at call sites.  A changed value writes + re-renders as before (user 2026-08-16:
//! providers self-report; every created grid announces, labeled, so stale never-used grids become visible).
void EmitAt(const std::string& section, const std::string& key, json value, Detail minLevel = Detail::Terse);

//! \brief PEAK resident set size of this process, in MB (Linux \c VmHWM: the high-water mark, not the
//! instantaneous \c VmRSS).  0 when unavailable.
//!
//! The HIGH-WATER mark is the number that matters for the CP2K comparison and for whether a cell fits the
//! box at all: an instantaneous reading taken after the streams are released, or between two peaks, tells
//! you nothing about what the run actually demanded (doc/OpenWork.md Step 1 — RAM is half that table, and
//! neither code was reporting it).
double PeakRSS_MB();

//! \brief Announce ONE symmetry FOLD, in the one format every fold site uses.
//!
//! WHY THIS EXISTS (doc/OpenWork.md Step 0b; standing user complaint): the folds are the largest
//! reductions in the code — 12–24× on a magnetic group, up to 48× cubic, 71× measured on the pair-stream
//! cache — and until now NONE of them printed anything.  A reduction nobody can see is a reduction nobody
//! trusts, and an ARMED fold is indistinguishable on the console from a disarmed one, which is exactly how
//! an opt-in 71× lever stays opt-in for months.  So every site that folds says so, in one line, ALWAYS —
//! including when it did NOT fold, because "no fold here" is the more actionable message.
//! \param site  where the fold acts ("XC mesh", "collocation streams", "V_loc {G}-star", …)
//! \param nOps  crystal ops in the fold (0 = not armed)
//! \param raw   items before folding; \a reps items actually evaluated.  Factor = raw/reps.
//! \param note  optional qualifier (e.g. "magnetic (Shubnikov)", "opt-in GPW_STREAM_FOLD").
void EmitFold(const std::string& site, size_t nOps, size_t raw, size_t reps,
              const std::string& note = std::string());

//! RAII: open a TOP-LEVEL section `name` (cursor descends into it) and RENDER it to
//! the console once, complete, when the scope closes (the incremental unit -- a crash
//! mid-assembly still renders what completed).  Use for a section ASSEMBLED across
//! code layers (basis); use EmitSection for one a single provider builds whole (scf).
class Section
{
public:
    explicit Section(const std::string& name);
    ~Section();
    Section(const Section&)            = delete;
    Section& operator=(const Section&) = delete;
private:
    std::string itsName;
    std::size_t itsDepth;
};

//! RAII: append a fresh row object to array `key` under the current cursor and descend
//! into it, so every Set(...) until this scope closes -- here AND in code this scope
//! calls into -- lands in that row.  Pops the cursor on close.
class Row
{
public:
    explicit Row(const std::string& key);
    ~Row();
    Row(const Row&)            = delete;
    Row& operator=(const Row&) = delete;
private:
    std::size_t itsDepth;
};

//============================================================================
// Renderers -- json -> representation.  The SINK (stdout, a file, a rolling log,
// an HDF5 group) is just the caller's ostream; the renderer stays GENERIC.
//============================================================================

//! Generic console renderer: layout INFERRED from json structure -- an object of
//! scalars renders as a two-column key/value table, an array of uniform objects as
//! a multi-column table (>=2 rows), anything more nested as an indented tree.  NO
//! per-section renderers; renders ANY subtree (that is what makes Emit incremental).
void RenderConsole(const json& node, std::ostream& os, Detail = Detail::Normal);

//! The trivial renderer: the model, verbatim (for the GUI/HDF5).
std::string RenderJson(const json& node);

//============================================================================
// Format hints -- field name -> unit + precision, consulted by RenderConsole so
// units live HERE (and in the schema doc), never baked into the data values.
// Unknown fields fall back to a %g-ish default.
//============================================================================
struct FormatHint { std::string unit; int precision; bool fixed = false; };  //!< fixed => decimal places, not sig-figs
FormatHint HintFor(const std::string& field);

} // export namespace qchem::report
