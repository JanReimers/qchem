// File: Common/Imp/RunControl.C  The process's write-once default control.
module;
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
module qchem.RunControl;

namespace qchem
{

namespace
{
//! THE PROCESS'S DEFAULT RunControl -- a mutable static, justified (CLAUDE.md "a mutable static needs a careful justification"):
//!  - WHAT: the command source a calculation polls when it was handed none.
//!  - WHY PER-PROCESS: what it serves is per-process -- a signal is delivered to a process, a runctl address is a pid, and a test binary
//!    is one process whose tests build their own facades.  User ruling 2026-11-25: "Allow, write-once": set ONCE by a test binary's
//!    main before any test runs (a second set throws, below); library code never sets it; rundeck and the GUI pass their own control.
//!  - THREAD SAFETY: written once under the mutex before any SCF runs, then only read (under the same mutex, once per iteration
//!    boundary: negligible).  The control's OWN Poll is the implementation's to make safe against its delivery side.
std::mutex                  g_processControlLock;
std::shared_ptr<RunControl> g_processControl;
bool                        g_processControlSet=false;
}

void SetProcessRunControl(std::shared_ptr<RunControl> control)
{
    std::lock_guard<std::mutex> g(g_processControlLock);
    if (g_processControlSet)
        throw std::logic_error("SetProcessRunControl: the process's default RunControl is WRITE-ONCE (user ruling 2026-11-25) and was already set; "
                               "a front end that needs its own control hands it to its calculations instead");
    g_processControl=std::move(control);
    g_processControlSet=true;
}

std::shared_ptr<RunControl> ProcessRunControl()
{
    std::lock_guard<std::mutex> g(g_processControlLock);
    return g_processControl;
}

std::string Name(RunCommand::Kind k)
{
    switch (k)
    {
    case RunCommand::Kind::Stop:        return "stop";
    case RunCommand::Kind::SaveAndStop: return "save-and-stop";
    case RunCommand::Kind::Checkpoint:  return "checkpoint";
    case RunCommand::Kind::SetThreads:  return "threads";
    }
    return "?";
}

} // namespace qchem
