// gtestmain.C Main entry point for all unit tests.

#include "gtest/gtest.h"
import qchem.SCFCap;     // an SCF that ends at its iteration cap FAILS the running test (user ruling 2026-10-29)
import qchem.FrontEnd.RunControl;   // run control for a test process: runctl, SIGTERM (a ctest TIMEOUT), the heartbeat (user 2026-11-25)
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <unistd.h>   // getpid (the default stopped-state name)

// There is no process-wide integrals cache (doc/IntegralCacheOwnershipPlan.md S6): every cache has an owner, so nothing to set up here.

namespace {
//! The heartbeat names the running test, so `scripts/runctl list` says which case a pid is on.
class CurrentTest : public testing::EmptyTestEventListener
{
public:
    explicit CurrentTest(qchem::frontend::FrontEndRunControl& c) : itsControl(c) {}
    void OnTestStart(const testing::TestInfo& t) override {itsControl.SetCurrent(std::string(t.test_suite_name())+"."+t.name());}
    void OnTestEnd  (const testing::TestInfo&)   override {itsControl.SetCurrent("");}
private:
    qchem::frontend::FrontEndRunControl& itsControl;
};
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    // EVERY SCF TEST FAILS IF AN SCF IN IT HITS NMaxIter UNCONVERGED (user 2026-10-29).  The SCF runs on the test's own
    // thread, so ADD_FAILURE lands in the running test.  A deliberately bounded loop says so with qchem::BoundedSCF.
    // The same hook fails a test whose SCF a RUN-CONTROL command stopped (qchem::ReportSCFStopped): never a silent pass.
    qchem::SetSCFCapHook([](const std::string& what){ ADD_FAILURE() << what; });

    // RUN CONTROL (user 2026-11-25: "for ctest runs it would be useful to be able to control ITMain"): this process's commands --
    // scripts/runctl <pid> ..., SIGTERM (a ctest TIMEOUT sends it: TIMEOUT_SIGNAL_NAME), Ctrl-C, SIGUSR1 -- reach every SCF through the
    // process's WRITE-ONCE default control (user ruling 2026-11-25, "Allow, write-once"), set here, ONCE, before any test runs.  The
    // control directory is the front end's (HOME is read here, never in the library).
    const char* h=std::getenv("HOME");
    const std::filesystem::path dir=std::filesystem::path(h ? h : ".")/"Code"/"qchem6-runs"/"control";
    auto control=std::make_shared<qchem::frontend::FrontEndRunControl>("ITMain", dir, (dir/("ITMain."+std::to_string(::getpid())+".stopped.h5")).string());
    control->InstallSignalHandlers();
    qchem::SetProcessRunControl(control);
    testing::UnitTest::GetInstance()->listeners().Append(new CurrentTest(*control));   // gtest owns the listener
    return RUN_ALL_TESTS();
}
