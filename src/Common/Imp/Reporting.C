// File: Common/Imp/Reporting.C  Implementation of qchem.Reporting.
//
//   Holds the global sink (a stack of run contexts + the keyed GlobalReport) and
//   the ONE generic console renderer.  The renderer infers layout from structure:
//   it walks a json subtree and, per node, chooses key/value table | multi-column
//   table | indented tree -- no section ever names its own layout.
module;
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>    // /proc/self/status (PeakRSS_MB)
#include <iostream>   // the [fold] announcements go to cout, not just the json record
#include <ostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <cassert>
#include <nlohmann/json.hpp>
#include "tabulate/table.hpp"

module qchem.Reporting;

using namespace tabulate;

namespace qchem::report {

//============================================================================
// Global sink state.  A stack of run contexts gives nesting for free: each frame
// carries the run's key and its (growing) json document; CurrentRunReport() is the
// top of the stack.  End() files the top frame into g_global under its key and pops.
//============================================================================
namespace {

// A cursor path segment: an object key OR an array index.  We store the PATH (not raw
// json* -- ordered_json is vector-backed, so an insert can reallocate and dangle a
// pointer) and re-resolve from the run root on each write; depth is tiny (<=3).
struct PathSeg
{
    bool        isIndex;
    std::string key;
    std::size_t idx;
    static PathSeg Key  (const std::string& k) { return { false, k, 0 }; }
    static PathSeg Index(std::size_t i)         { return { true, {}, i }; }
};

struct RunContext
{
    std::string          key;
    json                 doc = json::object();
    std::vector<PathSeg> savedCursor;            // the OUTER run's cursor, restored at End (nesting)
};

std::vector<RunContext> g_stack;                 // the run-nesting stack (empty == no run)
std::vector<PathSeg>    g_cursor;                // WHERE the next Set lands within the current run
json                    g_global = json::object();
std::ostream*           g_console = nullptr;      // set once at startup; null == silent
Detail                  g_detail  = Detail::Normal;

// Resolve the current cursor node within the current run's doc (null if no run open).
json* cursorNode()
{
    if (g_stack.empty()) return nullptr;
    json* node = &g_stack.back().doc;
    for (const auto& s : g_cursor) node = s.isIndex ? &(*node)[s.idx] : &(*node)[s.key];
    return node;
}

} // anonymous namespace

json& CurrentRunReport()
{
    assert(!g_stack.empty() && "CurrentRunReport() called outside a Begin/End bracket");
    return g_stack.back().doc;
}

const json& GlobalReport() { return g_global; }

void ClearGlobal()
{
    assert(g_stack.empty() && "ClearGlobal() while a run is open");
    g_global = json::object();
}

int Depth() { return static_cast<int>(g_stack.size()); }

bool InSection(const std::string& name)
{
    for (const auto& s : g_cursor) if (!s.isIndex && s.key == name) return true;
    return false;
}

void Begin(const std::string& name) { Begin(name, NowIso()); }

void Begin(const std::string& name, const std::string& startTime)
{
    RunContext ctx;
    ctx.key         = name + "@" + startTime;
    ctx.savedCursor = std::move(g_cursor);        // a nested run must not corrupt the outer cursor
    g_cursor.clear();                             // ...it starts at its own run root
    g_stack.push_back(std::move(ctx));
}

void End()
{
    assert(!g_stack.empty() && "End() with no matching Begin()");
    RunContext ctx = std::move(g_stack.back());
    g_stack.pop_back();
    g_cursor = std::move(ctx.savedCursor);        // restore the outer run's cursor
    g_global[ctx.key] = std::move(ctx.doc);       // every run kept (parent + nested + sequential)
}

std::string NowIso()
{
    // Production timestamp for Begin's startTime.  (Tests inject a fixed string, so
    // this wall-clock read never taints a golden-ish render.)
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream ss;
    ss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S");
    return ss.str();
}

void SetConsole(std::ostream& os, Detail d) { g_console = &os; g_detail = d; }
void ClearConsole()                          { g_console = nullptr; }

void Log(const std::string& message)
{
    if (g_console) *g_console << "· " << message << " …" << std::endl;   // endl => flush now (the heartbeat)
}

//============================================================================
// Scalar formatting -- a json leaf + its field's FormatHint -> a display string.
// Precision comes from the hint; the unit is appended for key/value leaves (table
// columns carry the unit in the header instead -- see renderTable).
//============================================================================
namespace {

bool isScalar(const json& v)                       // primitive leaf, or a flat list of them
{
    if (v.is_primitive()) return true;
    if (v.is_array())
        for (const auto& e : v) if (!e.is_primitive()) return false;
    return v.is_array();                            // empty array counts as a (trivial) scalar list
}

std::string fmtOne(const json& v, const FormatHint& h)
{
    if (v.is_null())            return "-";
    if (v.is_string())          return v.get<std::string>();
    if (v.is_boolean())         return v.get<bool>() ? "true" : "false";
    if (v.is_number_integer())  return std::to_string(v.get<long long>());
    if (v.is_number_unsigned()) return std::to_string(v.get<unsigned long long>());
    if (v.is_number_float())
    {
        std::ostringstream ss;
        if (h.fixed) ss << std::fixed;                     // decimal places (percentages, MB) vs sig-figs
        ss << std::setprecision(h.precision > 0 ? h.precision : 6) << v.get<double>();
        return ss.str();
    }
    return v.dump();
}

std::string fmtScalar(const std::string& field, const json& v, bool withUnit)
{
    const FormatHint h = HintFor(field);
    if (v.is_array())                               // flat list -> "a, b, c"
    {
        std::string s;
        for (std::size_t i = 0; i < v.size(); ++i) { if (i) s += ", "; s += fmtOne(v[i], h); }
        return s;
    }
    std::string s = fmtOne(v, h);
    if (withUnit && !h.unit.empty() && (v.is_number())) s += " " + h.unit;
    return s;
}

//============================================================================
// The generic renderer.  renderNode dispatches; renderKeyValue / renderTable do the
// two table shapes; deeper nesting recurses with indent (the "tree").
//============================================================================

std::string pad(int indent) { return std::string(std::size_t(indent) * 2, ' '); }

// A two-column key/value table (the degenerate depth-1 layout: an object of scalars).
void renderKeyValue(const std::vector<std::pair<std::string, const json*>>& rows,
                    std::ostream& os, int indent)
{
    Table t;
    for (const auto& [k, v] : rows)
        t.add_row(Table::Row_t{ k, fmtScalar(k, *v, /*withUnit*/ true) });
    t.format().hide_border();                       // key/value reads better without a grid
    std::istringstream in(t.str());
    for (std::string line; std::getline(in, line); ) os << pad(indent) << line << "\n";
}

// A multi-column table: an array of uniform objects (>=2 rows).  Columns are the
// union of keys in first-seen order; the unit (if any) rides in the header.
void renderTable(const json& arr, std::ostream& os, int indent)
{
    std::vector<std::string> cols;
    for (const auto& row : arr)
        for (auto it = row.begin(); it != row.end(); ++it)
            if (std::find(cols.begin(), cols.end(), it.key()) == cols.end()) cols.push_back(it.key());

    Table t;
    Table::Row_t header;
    for (const auto& c : cols)
    {
        const FormatHint h = HintFor(c);
        header.push_back(h.unit.empty() ? c : c + " (" + h.unit + ")");
    }
    t.add_row(header);

    for (const auto& row : arr)
    {
        Table::Row_t cells;
        for (const auto& c : cols)
            cells.push_back(row.contains(c) ? fmtScalar(c, row.at(c), /*withUnit*/ false) : "-");
        t.add_row(cells);
    }
    t.row(0).format().font_style({ FontStyle::bold });
    // Inter-row rules: if the LEADING column repeats down the table (a GROUPED table -- e.g. basis.usage,
    // many s rows then many p rows) keep a rule only BETWEEN groups (leading value changes) and none within,
    // so each group reads as one block.  If the leading column is all-distinct (e.g. perIrrep, one row per
    // irrep) drop every inter-row rule for a compact list.  Row 1's top (the header rule) is always kept.
    const std::string key0 = cols.empty() ? std::string() : cols.front();
    auto lead = [&](std::size_t r) -> std::string {                        // leading-column value of arr row r
        return arr[r].contains(key0) ? fmtScalar(key0, arr[r].at(key0), false) : std::string();
    };
    bool grouped = false;
    for (std::size_t r = 1; r < arr.size(); ++r) if (lead(r) == lead(r - 1)) { grouped = true; break; }
    for (std::size_t i = 2; i <= arr.size(); ++i)                          // tabulate row i == arr row i-1
    {
        const bool groupBoundary = grouped && lead(i - 1) != lead(i - 2);
        if (!groupBoundary) t.row(i).format().hide_border_top();           // rule only at a group boundary
    }

    std::istringstream in(t.str());
    for (std::string line; std::getline(in, line); ) os << pad(indent) << line << "\n";
}

// The field/section -> min console-detail map the RenderConsole skeleton anticipated: a few blocks are
// ALWAYS in the json but render only at higher verbosity (the raw basis-usage tuning data -- a consumer like
// valgen renders a friendlier joined view at Normal).  Unknown keys => Terse (shown at every level).
Detail keyMinDetail(const std::string& key)
{
    if (key == "exponents" || key == "usage") return Detail::Verbose;
    return Detail::Terse;
}

void renderNode(const json& node, std::ostream& os, int indent, Detail d)
{
    if (node.is_array())
    {
        // Array of >=2 uniform objects -> table; a single object (or empty) reads
        // better as a tree; a flat scalar list is one key/value-ish line.
        const bool objRows = !node.empty() &&
            std::all_of(node.begin(), node.end(), [](const json& e){ return e.is_object(); });
        if (objRows && node.size() >= 2) { renderTable(node, os, indent); return; }
        if (objRows)                                                     // single-row -> tree
        {
            for (const auto& e : node) renderNode(e, os, indent, d);
            return;
        }
        if (isScalar(node)) { os << pad(indent) << fmtScalar("", node, true) << "\n"; return; }
        for (const auto& e : node) renderNode(e, os, indent, d);         // ragged array -> recurse
        return;
    }

    if (node.is_object())
    {
        // Split scalar members (a key/value block) from container members (headed sub-nodes), skipping any
        // member gated above the current detail level.  Insertion order is preserved (ordered_json).
        std::vector<std::pair<std::string, const json*>> scalars;
        std::vector<std::pair<std::string, const json*>> containers;
        for (auto it = node.begin(); it != node.end(); ++it)
        {
            if ((int)keyMinDetail(it.key()) > (int)d) continue;          // Verbose-only block, below threshold
            (isScalar(it.value()) ? scalars : containers).emplace_back(it.key(), &it.value());
        }

        if (!scalars.empty()) renderKeyValue(scalars, os, indent);
        for (const auto& [k, v] : containers)
        {
            if (v->empty()) continue;                                    // skip empty sections
            os << pad(indent) << k << ":\n";
            renderNode(*v, os, indent + 1, d);
        }
        return;
    }

    os << pad(indent) << fmtScalar("", node, true) << "\n";              // bare scalar
}

} // anonymous namespace

void RenderConsole(const json& node, std::ostream& os, Detail d)
{
    renderNode(node, os, 0, d);                                          // d gates the field->min-level map
}

std::string RenderJson(const json& node) { return node.dump(2); }

// Render one completed section to the console -- incremental, depth-1 only (nested
// runs stay quiet on the console but are kept in full in the json).
static void renderSectionToConsole(const std::string& name, const json& section)
{
    if (g_console && Depth() == 1 && !section.empty())
    {
        *g_console << "\n" << name << "\n";
        RenderConsole(section, *g_console, g_detail);
    }
}

void EmitSection(const std::string& name, json section)
{
    CurrentRunReport()[name] = section;                    // ALWAYS record (json is complete)
    renderSectionToConsole(name, section);
}

void Set(const std::string& key, json value)
{
    if (json* n = cursorNode()) (*n)[key] = std::move(value);
}

void EmitAt(const std::string& section, const std::string& key, json value, Detail minLevel)
{
    if (g_stack.empty()) return;                       // inert outside a run
    // IDEMPOTENT: an identical value already at this path is a no-op (no rewrite, no re-render).  This is
    // what lets a provider announce itself UNCONDITIONALLY at its own trigger (e.g. every grid announces at
    // construction) with the dedup living HERE, run-scoped -- instead of in function-local-static latches at
    // the call sites, which are process-global and leak across runs (the old PWTerms void* latch).  A
    // CHANGED value still writes and re-renders, so late updates behave as before.
    {
        const json& doc = g_stack.back().doc;
        if (doc.contains(section) && doc[section].contains(key) && doc[section][key] == value) return;
    }
    g_stack.back().doc[section][key] = value;          // ALWAYS record (the json record is complete)
    // Console gate: render only at depth 1, non-empty, AND when the configured detail meets this block's
    // minimum (so a Verbose-only block like basis.usage stays in the json but prints only at Detail::Verbose).
    if (g_console && Depth() == 1 && !value.empty() && (int)g_detail >= (int)minLevel)
    {
        *g_console << "\n" << section << " ▸ " << key << "\n";
        RenderConsole(value, *g_console, g_detail);
    }
}

// ONE format for every fold site (see the header).  Goes to std::cout unconditionally -- a fold that only
// shows up in a json nobody opens is still invisible -- and into the run record under folds.<site> when a
// run is open.  Idempotent through EmitAt, so a site that re-announces the same fold prints once.
// VmHWM -- the kernel's own high-water mark, so it is correct no matter when we ask (unlike VmRSS, which
// reads whatever happens to be resident at the moment of the call).
double PeakRSS_MB()
{
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("VmHWM:", 0) == 0)
        {
            std::istringstream ss(line.substr(6));
            double kb = 0.0;
            ss >> kb;
            return kb / 1024.0;
        }
    return 0.0;
}

void EmitFold(const std::string& site, size_t nOps, size_t raw, size_t reps, const std::string& note)
{
    const double factor = reps ? double(raw)/double(reps) : 1.0;
    json j;
    j["ops"]=(long)nOps; j["raw"]=(long)raw; j["reps"]=(long)reps; j["factor"]=factor;
    if (!note.empty()) j["note"]=note;
    // RUN-SCOPED DEDUP, on the console as well as in the record.  A site announces UNCONDITIONALLY at its
    // own trigger -- but several triggers repeat per k-block (the V_loc {G}-star fires once per block), and
    // eight identical lines is noise that trains the reader to skip them.  Dedup lives here, run-scoped,
    // never in a process-global latch at the call site (see EmitAt's note).
    if (!g_stack.empty())
    {
        json& folds = g_stack.back().doc["folds"];
        if (folds.contains(site) && folds[site]==j) return;      // already announced, unchanged
        folds[site]=j;
    }
    // ARMED-NESS IS THE REDUCTION, not the op count: some sites know their fold factor but not how many
    // ops produced it (a Fold carries orbits, not the op set), so an absent nOps must not read as "off".
    std::cout << "[fold] " << site << ": ";
    if (reps>=raw)
        std::cout << "NONE (" << raw << " items evaluated in full)";
    else
    {
        if (nOps) std::cout << nOps << " ops, ";
        // std::defaultfloat restores the FORMAT flag but NOT the precision -- leaving cout at 2 for the rest
        // of the run, which silently truncated every energy printed after a fold line ("Etot=-7.1").  Restore
        // both, so a diagnostic can never degrade the numbers the run reports.
        const std::streamsize prec0=std::cout.precision();
        std::cout << raw << " -> " << reps << " representatives = "
                  << std::fixed << std::setprecision(2) << factor << "x" << std::defaultfloat;
        std::cout.precision(prec0);
    }
    if (!note.empty()) std::cout << "  [" << note << "]";
    std::cout << std::endl;
}

Section::Section(const std::string& name) : itsName(name), itsDepth(g_cursor.size())
{
    if (g_stack.empty()) return;                           // inert outside a run
    (*cursorNode())[name];                                 // ensure the member exists
    g_cursor.push_back(PathSeg::Key(name));
}

Section::~Section()
{
    g_cursor.resize(itsDepth);                             // pop back to the enclosing context
    if (g_stack.empty()) return;
    const json& doc = g_stack.back().doc;
    if (doc.contains(itsName)) renderSectionToConsole(itsName, doc[itsName]);
}

Row::Row(const std::string& key) : itsDepth(g_cursor.size())
{
    if (g_stack.empty()) return;                           // inert outside a run
    json& arr = (*cursorNode())[key];
    if (!arr.is_array()) arr = json::array();
    arr.push_back(json::object());
    g_cursor.push_back(PathSeg::Key(key));
    g_cursor.push_back(PathSeg::Index(arr.size() - 1));
}

Row::~Row() { g_cursor.resize(itsDepth); }

FormatHint HintFor(const std::string& field)
{
    // Field -> (unit, precision).  Small on purpose; grows with each migrated section.
    // Units: energies Hartree, lengths Bohr, exponents a.u.  Unknown -> {"", 6}.
    static const std::vector<std::pair<std::string, FormatHint>> hints = {
        { "densityEcut",   { "Ha", 6 } }, { "ecut",       { "Ha", 6 } },
        { "cond",          { "",   3 } }, { "lambdaMin",  { "",   3 } },
        { "lambdaMax",     { "",   3 } }, { "ramMB",      { "MB", 3, true } },
        { "reusePct",      { "%",  1, true } }, { "pct",   { "%",  1, true } },
        { "alpha",         { "a.u.", 4 } }, { "kappa",    { "",   4 } },
        { "cutoffFactor",  { "",   3 } }, { "fieldSharpness", { "", 4 } },
        { "smearingkT",    { "Ha", 6 } }, { "minDE",      { "",   2 } },
        { "minDrho",       { "",   2 } }, { "minDFD",     { "",   2 } },
        { "minFD",         { "",   2 } }, { "kerkerG0",   { "a.u.", 3 } },
        { "momSmearPenalty", { "Ha", 3 } },
    };
    for (const auto& [k, h] : hints) if (k == field) return h;
    return FormatHint{ "", 6 };
}

//============================================================================
// Timing ledger (see the interface doc).  Setup + the SCF loop are single-threaded (OMP parallelism
// lives INSIDE timed regions and never calls AddTime), so a plain map suffices -- same reasoning as
// the sink itself.  Insertion order is irrelevant; EmitTimings sorts by cost.
//============================================================================
namespace
{
std::map<std::string,double> g_times;
// CALL COUNTS beside the seconds.  A bucket's TOTAL is only comparable against another bucket's total
// when the two fire at the same CADENCE, and the ledger gave no way to know: "rho sampling 3.3 s" vs
// "matrix-free density 30.9 s" reads as a 9x only if both fire once per iteration, and they do not.
// (doc/OpenWork.md: the 1.70-vs-35.0 reading that made the low-rank win look 20x.)  One counter per
// bucket makes the per-call price -- the quantity an optimisation decision actually needs -- readable.
std::map<std::string,size_t> g_calls;
// Nesting stack for EXCLUSIVE times: each open Timed accumulates its children's ELAPSED here, and
// charges only (elapsed - children) to its own bucket -- so buckets are disjoint (see the interface doc).
// Setup + the SCF loop are single-threaded (OMP lives INSIDE timed regions), so a plain stack suffices.
std::vector<double> g_tchildren;
}

// COUNT ONLY THE ENTRIES THAT ACTUALLY DID WORK.  The time charged here is EXCLUSIVE (my elapsed minus my
// children's), so a scope whose body delegated to a nested Timed contributes ~0 seconds -- and counting it
// anyway divides real work by inflated entries and under-reports the per-call price.  MEASURED 2026-08-21:
// the XC rho-sampling bucket read 0.021 s/call over 83 entries when 42 of them had branched into the
// matrix-free child; the honest figure over the 41 entries that ran the GEMM was 0.042.  A bucket entered
// with nothing to charge is not a call, for the purpose of the question this number exists to answer.
void AddTime(const std::string& key, double seconds)
{
    g_times[key]+=seconds;
    if (seconds > 1e-9) g_calls[key]++;    // sub-ns == a pass-through entry, not a unit of work
}

Timed::Timed(std::string key)
    : itsKey(std::move(key))
    , itsT0(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now().time_since_epoch()).count())
{
    g_tchildren.push_back(0.0);
}

Timed::~Timed()
{
    const long long t1=std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now().time_since_epoch()).count();
    const double elapsed=double(t1-itsT0)*1e-9;
    AddTime(itsKey, elapsed-g_tchildren.back());          // exclusive: my elapsed minus my children's
    g_tchildren.pop_back();
    if (!g_tchildren.empty()) g_tchildren.back()+=elapsed; // I am my parent's child, whole
}

void EmitTimings(const std::string& name)
{
    if (g_times.empty()) return;
    std::vector<std::pair<std::string,double>> rows(g_times.begin(), g_times.end());
    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b){return a.second>b.second;});
    json j;
    // The per-call price rides in the LABEL, not as a nested object: a sub-object per bucket would turn
    // the flat sorted-by-cost table into a tree (renderNode heads any container member), and the flat
    // table is the whole readability of this report.  Suffix only when the bucket fired more than once --
    // a one-shot setup bucket has no cadence to report and the noise would cost more than it tells.
    for (const auto& [k,s] : rows)
    {
        const size_t n=g_calls.count(k) ? g_calls.at(k) : 0;
        if (n>1)
        {
            std::ostringstream lbl;
            lbl<<k<<"  [x"<<n<<", "<<std::setprecision(4)<<(s/double(n))<<" s/call]";
            j[lbl.str()]=s;
        }
        else j[k]=s;                            // ordered_json keeps the sorted-by-cost order
    }
    g_times.clear();
    g_calls.clear();
    // PEAK RAM beside the time buckets: the two halves of "what did this run cost" belong in one place,
    // and until now the report answered only one of them (doc/OpenWork.md Step 1).  VmHWM is the
    // process-wide high-water mark, so on a multi-run process it is monotone -- it reports the
    // WATERMARK SO FAR, which is the honest reading for a benchmark row run one config per process.
    if (const double mb=PeakRSS_MB(); mb>0.0) j["PEAK RSS (MB, process high-water)"]=mb;
    EmitSection(name, std::move(j));
}

} // namespace qchem::report
