// File: Common/Imp/Reporting.C  Implementation of qchem.Reporting.
//
//   Holds the global sink (a stack of run contexts + the keyed GlobalReport) and
//   the ONE generic console renderer.  The renderer infers layout from structure:
//   it walks a json subtree and, per node, chooses key/value table | multi-column
//   table | indented tree -- no section ever names its own layout.
module;
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
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

    std::istringstream in(t.str());
    for (std::string line; std::getline(in, line); ) os << pad(indent) << line << "\n";
}

void renderNode(const json& node, std::ostream& os, int indent)
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
            for (const auto& e : node) renderNode(e, os, indent);
            return;
        }
        if (isScalar(node)) { os << pad(indent) << fmtScalar("", node, true) << "\n"; return; }
        for (const auto& e : node) renderNode(e, os, indent);            // ragged array -> recurse
        return;
    }

    if (node.is_object())
    {
        // Split scalar members (a key/value block) from container members (headed
        // sub-nodes).  Insertion order is preserved (ordered_json).
        std::vector<std::pair<std::string, const json*>> scalars;
        std::vector<std::pair<std::string, const json*>> containers;
        for (auto it = node.begin(); it != node.end(); ++it)
            (isScalar(it.value()) ? scalars : containers).emplace_back(it.key(), &it.value());

        if (!scalars.empty()) renderKeyValue(scalars, os, indent);
        for (const auto& [k, v] : containers)
        {
            if (v->empty()) continue;                                    // skip empty sections
            os << pad(indent) << k << ":\n";
            renderNode(*v, os, indent + 1);
        }
        return;
    }

    os << pad(indent) << fmtScalar("", node, true) << "\n";              // bare scalar
}

} // anonymous namespace

void RenderConsole(const json& node, std::ostream& os, Detail /*d*/)
{
    // Detail-level filtering is a later step (a field/section -> min-level map on this
    // side); the skeleton renders the whole subtree.
    renderNode(node, os, 0);
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

} // namespace qchem::report
