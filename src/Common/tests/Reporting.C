// File: src/Common/tests/Reporting.C  Tests for qchem.Reporting (the schema-check +
//   render tests that stand in for the compile-time key-safety we trade for json).
#include "gtest/gtest.h"
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>
import qchem.Reporting;

using namespace qchem::report;

// Hand-build a report shaped like the doc/RunReportPlan.md schema and assert its
// shape.  This IS the "compile-time-check replacement": if a writer helper later
// produces the wrong key/shape, a test like this catches it.
static json SampleRun()
{
    json r = json::object();
    r["meta"]  = { { "title", "H2O / LDA" }, { "nElectrons", 8 }, { "spinPolarized", false } };
    r["scf"]   = { { "standard", { { "nMaxIter", 50 }, { "mixer", "Ker" }, { "smearingkT", 0.0 } } } };
    r["basis"] = { { "name", "valence_lowq" }, { "nFunctions", 23 },
                   { "perIrrep", json::array({
                        { { "irrep", "A1" }, { "nFunctions", 9 }, { "cond", 12.3 } },
                        { { "irrep", "B2" }, { "nFunctions", 5 }, { "cond", 4.5 } } }) } };
    return r;
}

// -- schema shape -----------------------------------------------------------

TEST(Reporting, SchemaShape)
{
    const json r = SampleRun();
    ASSERT_TRUE(r.contains("meta"));
    ASSERT_TRUE(r.contains("basis"));
    EXPECT_EQ(r["meta"]["nElectrons"].get<int>(), 8);
    ASSERT_TRUE(r["basis"]["perIrrep"].is_array());
    EXPECT_EQ(r["basis"]["perIrrep"].size(), 2u);
    EXPECT_EQ(r["basis"]["perIrrep"][0]["irrep"].get<std::string>(), "A1");
}

// -- ordered_json keeps insertion order (sections must not sort alphabetically) --

TEST(Reporting, SectionsKeepInsertionOrder)
{
    json r = json::object();
    r["zeta"] = 1; r["alpha"] = 2; r["mu"] = 3;
    std::vector<std::string> keys;
    for (auto it = r.begin(); it != r.end(); ++it) keys.push_back(it.key());
    EXPECT_EQ(keys, (std::vector<std::string>{ "zeta", "alpha", "mu" }));
}

// -- generic renderer: layout inferred from structure -----------------------

TEST(Reporting, RenderKeyValueBlock)
{
    json meta = { { "title", "H2O" }, { "nElectrons", 8 } };
    std::ostringstream os;
    RenderConsole(meta, os);
    const std::string out = os.str();
    EXPECT_NE(out.find("title"), std::string::npos);
    EXPECT_NE(out.find("H2O"), std::string::npos);
    EXPECT_NE(out.find("nElectrons"), std::string::npos);
    EXPECT_NE(out.find('8'), std::string::npos);
}

TEST(Reporting, RenderArrayOfObjectsAsTable)
{
    // >=2 uniform objects -> a multi-column table: every irrep + column appears.
    json perIrrep = json::array({
        { { "irrep", "A1" }, { "nFunctions", 9 } },
        { { "irrep", "B2" }, { "nFunctions", 5 } } });
    std::ostringstream os;
    RenderConsole(perIrrep, os);
    const std::string out = os.str();
    EXPECT_NE(out.find("irrep"), std::string::npos);
    EXPECT_NE(out.find("nFunctions"), std::string::npos);
    EXPECT_NE(out.find("A1"), std::string::npos);
    EXPECT_NE(out.find("B2"), std::string::npos);
}

TEST(Reporting, EmptySectionsAreSkipped)
{
    json r = { { "populated", { { "k", "v" } } }, { "empty", json::object() } };
    std::ostringstream os;
    RenderConsole(r, os);
    EXPECT_EQ(os.str().find("empty"), std::string::npos);
    EXPECT_NE(os.str().find("populated"), std::string::npos);
}

// -- the sink: Begin/End keying, nesting, Depth -----------------------------

TEST(Reporting, BeginEndKeysAndFilesTheRun)
{
    ClearConsole();
    EXPECT_EQ(Depth(), 0);
    Begin("H2O", "2026-07-26T10:00:00");
    EXPECT_EQ(Depth(), 1);
    EmitSection("meta", { { "title", "H2O" } });
    EXPECT_TRUE(CurrentRunReport().contains("meta"));
    End();
    EXPECT_EQ(Depth(), 0);
    EXPECT_TRUE(GlobalReport().contains("H2O@2026-07-26T10:00:00"));
    EXPECT_EQ(GlobalReport()["H2O@2026-07-26T10:00:00"]["meta"]["title"].get<std::string>(), "H2O");
}

TEST(Reporting, BeginAutoStampsTimestamp)
{
    ClearConsole();
    Begin("Auto");                                  // no timestamp supplied
    EmitSection("meta", { { "k", "v" } });
    End();
    // Exactly one key "Auto@<something-nonempty>" exists.
    bool found = false;
    for (auto it = GlobalReport().begin(); it != GlobalReport().end(); ++it)
        if (it.key().rfind("Auto@", 0) == 0 && it.key().size() > std::string("Auto@").size())
            found = true;
    EXPECT_TRUE(found);
}

TEST(Reporting, NestedRunsDoNotClobber)
{
    ClearConsole();
    Begin("outer", "t1");
    EmitSection("meta", { { "which", "outer" } });
    EXPECT_EQ(Depth(), 1);

    Begin("innerSAD", "t2");                       // the HF/DHF SAD-bootstrap shape
    EXPECT_EQ(Depth(), 2);
    EmitSection("meta", { { "which", "inner" } });
    EXPECT_EQ(CurrentRunReport()["meta"]["which"].get<std::string>(), "inner");
    End();

    EXPECT_EQ(Depth(), 1);                          // back in the outer run, unclobbered
    EXPECT_EQ(CurrentRunReport()["meta"]["which"].get<std::string>(), "outer");
    End();

    EXPECT_TRUE(GlobalReport().contains("outer@t1"));
    EXPECT_TRUE(GlobalReport().contains("innerSAD@t2"));   // both runs kept
}

TEST(Reporting, ConsoleRendersOnlyAtDepth1)
{
    std::ostringstream os;
    SetConsole(os, Detail::Normal);
    Begin("outer", "t1");
    EmitSection("a", { { "x", 1 } });               // depth 1 -> renders
    Begin("inner", "t2");
    EmitSection("b", { { "y", 2 } });               // depth 2 -> silent on console...
    End();
    End();
    ClearConsole();
    EXPECT_NE(os.str().find("a"), std::string::npos);
    EXPECT_EQ(os.str().find('y'), std::string::npos);   // ...but 'b'/y never hit the console
    // both are still in the json record:
    EXPECT_TRUE(GlobalReport()["inner@t2"].contains("b"));
}

// -- a provider-built section round-trips through EmitSection ----------------
// (The provider owns the keys -- here we stand in for one: build the section json,
// emit it, and confirm the sink stored exactly what was built.  The REAL scf-schema
// check lives with its provider, in IntegrationTests/M_Calculation.C.)

TEST(Reporting, ProviderBuildsSectionJsonDirectly)
{
    ClearConsole();
    Begin("scfTest", "t");
    json scf;                                     // a provider building json by hand
    scf["standard"]["nMaxIter"] = 50;
    scf["standard"]["mixer"]    = "Ker";
    scf["advanced"]["guard"]["maxReleases"] = 2;
    EmitSection("scf", scf);
    const json& stored = CurrentRunReport()["scf"];
    End();

    EXPECT_EQ(stored, scf);                       // stored verbatim, no struct in between
    EXPECT_EQ(stored["standard"]["nMaxIter"].get<long>(), 50);
    EXPECT_EQ(stored["advanced"]["guard"]["maxReleases"].get<int>(), 2);
}

// -- the section cursor (context-free deep writes) --------------------------

TEST(Reporting, CursorSetLandsAtSectionAndRow)
{
    ClearConsole();
    Begin("cur", "t");
    {
        Section basis("basis");
        Set("name", "dzvp");                     // lands at basis level
        {
            Row r("perIrrep");
            Set("irrep", "A1");
            Set("cond", 12.3);                   // stand-in for a deep provider's write
        }
        {
            Row r("perIrrep");
            Set("irrep", "B2");
            Set("cond", 4.5);
        }
    }
    const json& b = CurrentRunReport()["basis"];
    End();

    EXPECT_EQ(b["name"].get<std::string>(), "dzvp");
    ASSERT_TRUE(b["perIrrep"].is_array());
    ASSERT_EQ(b["perIrrep"].size(), 2u);
    EXPECT_EQ(b["perIrrep"][0]["irrep"].get<std::string>(), "A1");
    EXPECT_EQ(b["perIrrep"][0]["cond"].get<double>(), 12.3);
    EXPECT_EQ(b["perIrrep"][1]["irrep"].get<std::string>(), "B2");
}

TEST(Reporting, CursorInertWithoutRun)
{
    ClearConsole();
    // No Begin: every cursor op must be a safe no-op (kernels used outside a report).
    Section s("basis");
    Set("name", "dzvp");
    Row r("perIrrep");
    Set("cond", 1.0);
    SUCCEED();                                    // reaching here without an assert/crash IS the test
}

TEST(Reporting, CursorRestoredAfterScopeAndException)
{
    ClearConsole();
    Begin("cur", "t");
    try
    {
        Section basis("basis");
        Row r("perIrrep");
        Set("irrep", "A1");
        throw std::runtime_error("mid-assembly boom");   // unwinds through Row + Section
    }
    catch (...) {}
    // Cursor unwound cleanly: a fresh top-level section is still addressable.
    { Section scf("scf"); Set("nMaxIter", 20); }
    const json doc = CurrentRunReport();          // COPY before End (End moves the doc out)
    End();
    EXPECT_EQ(doc["basis"]["perIrrep"][0]["irrep"].get<std::string>(), "A1");  // partial kept
    EXPECT_EQ(doc["scf"]["nMaxIter"].get<int>(), 20);                          // cursor not corrupted
}

TEST(Reporting, NestedRunCursorIsolated)
{
    ClearConsole();
    Begin("outer", "t1");
    Section basis("basis");
    Set("name", "outerBasis");                   // outer cursor at basis
    {
        Begin("inner", "t2");                    // nested run: outer cursor must be preserved
        { Section g("grids"); Set("raster", 5); }
        End();
    }
    Set("engine", "mnd");                        // still lands in outer basis, cursor intact
    const json outer = CurrentRunReport();       // COPY before End (End moves the doc out)
    End();
    EXPECT_EQ(outer["basis"]["name"].get<std::string>(), "outerBasis");
    EXPECT_EQ(outer["basis"]["engine"].get<std::string>(), "mnd");
    EXPECT_TRUE(GlobalReport()["inner@t2"]["grids"].contains("raster"));
}

// -- RenderJson is the model verbatim ---------------------------------------

TEST(Reporting, RenderJsonRoundTrips)
{
    const json r = SampleRun();
    const std::string dumped = RenderJson(r);
    EXPECT_EQ(json::parse(dumped), r);
}

// Visual eyeball of the console layout (run: --gtest_also_run_disabled_tests
// --gtest_filter='*VisualDump*').  Not an assertion -- a human-readable sample.
TEST(Reporting, DISABLED_VisualDump)
{
    std::cout << "\n===== whole-run RenderConsole =====\n";
    RenderConsole(SampleRun(), std::cout);
    std::cout << "\n===== incremental EmitSection (console) =====\n";
    SetConsole(std::cout, Detail::Normal);
    Begin("H2O", "2026-07-26T10:00:00");
    const json run = SampleRun();
    for (auto it = run.begin(); it != run.end(); ++it) EmitSection(it.key(), it.value());
    End();
    ClearConsole();
}
