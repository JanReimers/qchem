// File: IntegrationTests/RunControlFrontEnd.C  The run-control FRONT END ITMain and rundeck share (qchem.FrontEnd.RunControl): the
// control-file grammar, the file channel (take, apply, log), the heartbeat, and the sticky stop.  No SCF: the facade half of the claim is
// UTCalculation's FacadeMatrix cell RunControlStopsAndSaves.  (OpenWork "Run control", user 2026-11-25.)
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
import qchem.FrontEnd.RunControl;

using namespace qchem;
using qchem::frontend::ParseControlLine;
using qchem::frontend::FrontEndRunControl;
namespace fs=std::filesystem;

TEST(FrontEndRunControl, EveryControlLineParsesOrSaysWhy)
{
    RunCommand c; std::string why;
    ASSERT_TRUE(ParseControlLine("threads 8", c, why)) << why;
    EXPECT_EQ(c.kind, RunCommand::Kind::SetThreads); EXPECT_EQ(c.threads.openmp, 8); EXPECT_EQ(c.threads.blas, 0);
    ASSERT_TRUE(ParseControlLine("  threads 4 blas 2 ", c, why)) << why;
    EXPECT_EQ(c.threads.openmp, 4); EXPECT_EQ(c.threads.blas, 2);
    ASSERT_TRUE(ParseControlLine("save-and-stop", c, why)) << why;
    EXPECT_EQ(c.kind, RunCommand::Kind::SaveAndStop); EXPECT_TRUE(c.path.empty());
    ASSERT_TRUE(ParseControlLine("checkpoint /tmp/x.h5", c, why)) << why;
    EXPECT_EQ(c.kind, RunCommand::Kind::Checkpoint); EXPECT_EQ(c.path, "/tmp/x.h5");
    ASSERT_TRUE(ParseControlLine("stop", c, why)) << why;
    EXPECT_EQ(c.kind, RunCommand::Kind::Stop);
    for (const char* bad : {"threads", "threads 0", "threads 4 openmp 2", "threads 4 blas 0", "halt", "stop now"})
    {
        why.clear();
        EXPECT_FALSE(ParseControlLine(bad, c, why)) << "'" << bad << "' parsed";
        EXPECT_FALSE(why.empty()) << "'" << bad << "' refused without saying why";
    }
}

TEST(FrontEndRunControl, TheControlFileIsTakenAppliedLoggedAndAStopIsSticky)
{
    const fs::path dir=fs::path(testing::TempDir())/("runctl_"+std::to_string(::getpid()));
    fs::remove_all(dir);
    {
        FrontEndRunControl control("test", dir, (dir/"default.h5").string());
        const fs::path progress=dir/(std::to_string(::getpid())+".progress");
        ASSERT_TRUE(fs::exists(progress)) << "no heartbeat at construction";
        EXPECT_TRUE(control.Poll().empty()) << "no file, no signal: nothing";
        { std::ofstream o(control.ControlFile()); o<<"threads 3\n# a comment\nbogus line\ncheckpoint\n"; }
        auto cmds=control.Poll();
        ASSERT_EQ(cmds.size(), 2u) << "threads + checkpoint (the comment skipped, the bad line ignored)";
        EXPECT_EQ(cmds[0].kind, RunCommand::Kind::SetThreads);
        EXPECT_EQ(cmds[1].kind, RunCommand::Kind::Checkpoint);
        EXPECT_FALSE(fs::exists(control.ControlFile())) << "the control file was not taken";
        std::ifstream done(control.ControlFile().string()+".done"); std::stringstream log; log<<done.rdbuf();
        EXPECT_NE(log.str().find("threads 3   # applied"), std::string::npos) << log.str();
        EXPECT_NE(log.str().find("bogus line   # IGNORED"), std::string::npos) << log.str();
        EXPECT_EQ(control.DefaultSavePath(), (dir/"default.h5").string());
        { std::ofstream o(control.ControlFile()); o<<"stop\n"; }
        cmds=control.Poll();
        ASSERT_EQ(cmds.size(), 1u); EXPECT_EQ(cmds[0].kind, RunCommand::Kind::Stop);
        cmds=control.Poll();
        ASSERT_EQ(cmds.size(), 1u) << "a stop is STICKY: every later poll answers stop";
        EXPECT_EQ(cmds[0].kind, RunCommand::Kind::Stop);
        std::ifstream hb(progress); std::stringstream s; s<<hb.rdbuf();
        EXPECT_NE(s.str().find("\"stopped\": true"), std::string::npos) << s.str();
        EXPECT_NE(s.str().find("\"polls\": 4"), std::string::npos) << s.str();
    }
    EXPECT_FALSE(fs::exists(dir/(std::to_string(::getpid())+".progress"))) << "a clean exit leaves its heartbeat behind";
    fs::remove_all(dir);
}
