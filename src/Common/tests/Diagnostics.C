// File: src/Common/tests/Diagnostics.C  qchem.Diagnostics (D-ENV step 2): the registry, the parse, the typo catch.
#include "gtest/gtest.h"
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

import qchem.Diagnostics;
namespace D=qchem::Diagnostics;

TEST(Diagnostics, ParseListTakesIdsAndValuesAndReportsTypos)
{
    std::vector<std::string> unknown;
    const auto m=D::ParseList("dm_rank,mesh_ortho=4,rss_trce,list", unknown);
    EXPECT_EQ(m.at("dm_rank"), "1") << "a bare id is switched on with the value 1";
    EXPECT_EQ(m.at("mesh_ortho"), "4") << "id=value carries the argument";
    EXPECT_TRUE(m.count("list"));
    ASSERT_EQ(unknown.size(), 1u);
    EXPECT_EQ(unknown[0], "rss_trce") << "a mistyped id is reported, not silently ignored";
    EXPECT_FALSE(m.count("rss_trce"));
}

TEST(Diagnostics, EveryRegisteredIdIsUniqueAndDescribed)
{
    std::set<std::string> ids, envs;
    for (const auto& e : D::Registry())
    {
        EXPECT_TRUE(ids.insert(e.id).second) << "duplicate id " << e.id;
        if (!e.legacyEnv.empty()) EXPECT_TRUE(envs.insert(e.legacyEnv).second) << "duplicate legacy name " << e.legacyEnv;
        EXPECT_FALSE(e.what.empty()) << e.id << " has no description";
    }
    std::ostringstream os; D::Describe(os);
    EXPECT_NE(os.str().find("dm_rank"), std::string::npos);
}

TEST(Diagnostics, AnUnregisteredIdInSourceIsABugNotAQuietOff)
{
    EXPECT_THROW(D::Enabled("no_such_diagnostic"), std::logic_error);
    EXPECT_THROW(D::Scoped("no_such_diagnostic"), std::logic_error);
}

TEST(Diagnostics, ScopedSwitchesOneIdForATestAndRestores)
{
    ASSERT_FALSE(D::Enabled("rho_negative")) << "the suite runs with diagnostics off";
    {
        D::Scoped on("rho_negative");
        EXPECT_TRUE(D::Enabled("rho_negative"));
        EXPECT_EQ(*D::Value("rho_negative"), "1");
        D::Scoped arg("mesh_ortho", "4");
        EXPECT_EQ(*D::Value("mesh_ortho"), "4");
    }
    EXPECT_FALSE(D::Enabled("rho_negative"));
    EXPECT_FALSE(D::Value("mesh_ortho").has_value());
}
