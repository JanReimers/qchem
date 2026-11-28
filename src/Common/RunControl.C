// File: Common/RunControl.C  Commands to a RUNNING calculation: stop, save-and-stop, checkpoint, change the thread count.
//
// OpenWork "Run control: commands to a running calculation" (user 2026-11-01, re-asked 2026-11-25: "it would still be nice to have a way
// to send signals to rundeck: 1) change thread count, 2) stop and save after the next iter completes").  Design: the lane G note
// laneG-runcontrol-design.md, approved by the controller 2026-11-25.
//
// THE SPLIT.  The library sees only this abstract face, and only POLLS it, at ONE safe point: the end of an SCF iteration (never mid
// Fock build).  DELIVERY is the front end's business: rundeck maps signals and a control file onto it, ITMain's main does the same for
// a test process, the GUI bridge (later) pushes into a queue.  What a command MEANS (save where, stop which stages) is the policy layer's:
// tCalculation acts on it (qchem.Calculation).
module;
#include <memory>
#include <string>
#include <vector>
export module qchem.RunControl;
export import qchem.Parallel;   // ThreadCounts: what a SetThreads carries

export namespace qchem
{

//! \brief One command to a running calculation.  Applied at the next iteration boundary, in the order received.
struct RunCommand
{
    enum class Kind
    {
        Stop,          //!< end the SCF after the iteration in progress; save nothing beyond what the run saves anyway
        SaveAndStop,   //!< save the state (\c path, else the run's \c saveStateTo), then stop as \c Stop
        Checkpoint,    //!< save the state now (\c path, else \c saveStateTo) and carry on
        SetThreads     //!< from the next iteration on, run on \c threads (OpenMP and the BLAS pool; 0 keeps that pool's count)
    };
    Kind         kind = Kind::Stop;
    ThreadCounts threads{};    //!< SetThreads only
    std::string  path;         //!< SaveAndStop / Checkpoint: where to save; empty = the run's own \c saveStateTo, else the control's default
};

//! \brief THE FACE a front end implements and a calculation polls.  \c Poll drains the commands that arrived since the last call; it
//! is called ONLY at an SCF iteration boundary, on the SCF's own thread, so an implementation needs to be safe only against its OWN
//! delivery side (a signal handler, a GUI thread).
class RunControl
{
public:
    virtual ~RunControl() = default;
    virtual std::vector<RunCommand> Poll() = 0;
    //! Where a save-and-stop or checkpoint with no \c path goes when the run has no \c saveStateTo either: the front end's default
    //! (rundeck: its run directory; a test binary: the control directory).  Empty = nowhere, and the stop says no state was saved.
    virtual std::string DefaultSavePath() const {return "";}
};

//! \brief THE PROCESS'S DEFAULT control: what a calculation that was handed none polls.  WRITE-ONCE (user ruling 2026-11-25, "Allow,
//! write-once"): a test binary's main sets it ONCE, before any test runs, so the hundreds of tests that build their own facades are
//! controllable (runctl, the ctest TIMEOUT signal) without each passing a control.  A SECOND set throws; library code never sets it;
//! rundeck and the GUI always hand their calculations their OWN control and never rely on it.  The mutable static behind it is
//! justified at its definition (Imp/RunControl.C).
void SetProcessRunControl(std::shared_ptr<RunControl> control);
//! The process's default control, or null when none was set (every production driver).
std::shared_ptr<RunControl> ProcessRunControl();

//! The name a report or a console line gives a command kind ("stop", "save-and-stop", "checkpoint", "threads").
std::string Name(RunCommand::Kind k);

} // namespace qchem
