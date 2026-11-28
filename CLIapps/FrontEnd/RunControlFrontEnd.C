// File: CLIapps/FrontEnd/RunControlFrontEnd.C  The run-control FRONT END shared by rundeck and ITMain (qchem.RunControl's deliverer).
//
// OpenWork "Run control" (user 2026-11-25); design laneG-runcontrol-design.md, approved 2026-11-25.  Everything here is the front end's
// business and never the library's: POSIX signals, the control directory under $HOME, the pid address, the heartbeat file.  The library
// sees only qchem::RunControl::Poll at an SCF iteration boundary.
//
// THE ADDRESS: one directory, ~/Code/qchem6-runs/control/, keyed by pid (an ITMain test has no run directory):
//   <pid>.progress   the HEARTBEAT, rewritten at every poll (write-temp-then-rename, so a reader never sees half a file)
//   <pid>.control    the COMMANDS, one per line, written by scripts/runctl; consumed at the next poll and moved to <pid>.control.done
// THE SIGNALS: SIGTERM and SIGINT (and a ctest TIMEOUT's TIMEOUT_SIGNAL_NAME SIGTERM) = save-and-stop, and the SECOND one kills (the
// handler is SA_RESETHAND); SIGUSR1 = save-and-stop, repeatable.  A stop is STICKY for the process: every later poll answers stop, so
// a test binary's later SCFs end at their first boundary and a deck's later stages never start.
module;
#include <atomic>
#include <chrono>
#include <csignal>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>   // getpid
export module qchem.FrontEnd.RunControl;
export import qchem.RunControl;

export namespace qchem::frontend
{

//! \brief Parse ONE control-file line into a command: `threads N [blas M]`, `save-and-stop [PATH]`, `checkpoint [PATH]`, `stop`.
//! FALSE (and \a why says what is wrong) for anything else; a bad line is reported and skipped, never fatal.
bool ParseControlLine(const std::string& line, RunCommand& out, std::string& why);

//! \brief rundeck's and ITMain's RunControl: the signal flag, the pid-keyed control file, the heartbeat.  ONE per process (the signal
//! handlers are per process); construct it in \c main.
class FrontEndRunControl : public RunControl
{
public:
    //! \a who: what this process is, for `runctl list` ("rundeck NaF.json -> runs/NaF", "ITMain"); \a dir: the control directory;
    //! \a defaultSave: where a save with no other path goes (rundeck: its run dir's states/; ITMain: \a dir).
    FrontEndRunControl(std::string who, std::filesystem::path dir, std::string defaultSave);
    ~FrontEndRunControl() override;   //!< removes this pid's heartbeat (a clean exit leaves nothing for runctl list to prune)
    FrontEndRunControl(const FrontEndRunControl&)=delete;
    FrontEndRunControl& operator=(const FrontEndRunControl&)=delete;

    //! Install the SIGTERM / SIGINT / SIGUSR1 handlers (only the atomic flag below is touched in them).
    void InstallSignalHandlers();
    //! What the process is working on now (ITMain: the running test's name); shown by the heartbeat.
    void SetCurrent(const std::string& what);

    std::vector<RunCommand> Poll() override;
    std::string DefaultSavePath() const override {return itsDefaultSave;}

    const std::filesystem::path& ControlFile() const {return itsControl;}
private:
    void WriteHeartbeat(const std::string& lastCommand);
    std::string           itsWho;
    std::filesystem::path itsDir, itsControl, itsProgress;
    std::string           itsDefaultSave;
    std::mutex            itsLock;          //!< SetCurrent (a gtest listener) vs Poll (the SCF thread)
    std::string           itsCurrent;
    size_t                itsPolls=0;
    bool                  itsStopped=false; //!< sticky: once a stop went out, every poll answers stop
    std::time_t           itsStarted;
};

} // namespace qchem::frontend

module :private;

namespace qchem::frontend
{

namespace
{
//! THE SIGNAL FLAG -- the one mutable static a signal handler may touch (CLAUDE.md, mutable statics; design note (a), approved 2026-11-25):
//!  - WHAT: a stop request a signal delivered (0 none, 1 save-and-stop) not yet seen by Poll.
//!  - WHY PER-PROCESS: a POSIX signal is delivered to the process and its handler can reach no object; async-signal-safety allows only a
//!    lock-free atomic store (a mutex, an allocation or a stream in a handler is undefined behaviour).
//!  - THREAD SAFETY: a lock-free std::atomic, set by the handler and exchanged to 0 by FrontEndRunControl::Poll on the SCF thread.
//!  It lives in the front end only (rundeck, ITMain); the library never sees it.
std::atomic<int> g_signalRequest{0};
static_assert(std::atomic<int>::is_always_lock_free, "a signal handler may only touch a lock-free atomic");

void OnStopSignal(int) {g_signalRequest.store(1);}   // async-signal-safe: one lock-free atomic store, nothing else

std::string Trim(std::string s)
{
    const auto b=s.find_first_not_of(" \t\r\n");
    if (b==std::string::npos) return "";
    const auto e=s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

std::string JsonString(const std::string& s)
{
    std::string o="\"";
    for (char c : s) { if (c=='"' || c=='\\') o+='\\'; if (c=='\n') { o+="\\n"; continue; } o+=c; }
    return o+"\"";
}
} // namespace

bool ParseControlLine(const std::string& raw, RunCommand& out, std::string& why)
{
    std::istringstream is(Trim(raw));
    std::string verb; is>>verb;
    out=RunCommand{};
    if (verb=="stop")          { out.kind=RunCommand::Kind::Stop; }
    else if (verb=="save-and-stop") { out.kind=RunCommand::Kind::SaveAndStop; is>>out.path; }
    else if (verb=="checkpoint")    { out.kind=RunCommand::Kind::Checkpoint;  is>>out.path; }
    else if (verb=="threads")
    {
        out.kind=RunCommand::Kind::SetThreads;
        int n=0; if (!(is>>n) || n<1) { why="threads needs a positive OpenMP count: `threads N [blas M]`"; return false; }
        out.threads.openmp=n;
        std::string b;
        if (is>>b)
        {
            int m=0;
            if (b!="blas" || !(is>>m) || m<1) { why="after `threads N` only `blas M` (M >= 1) may follow"; return false; }
            out.threads.blas=m;
        }
    }
    else { why="unknown command '"+verb+"' (stop | save-and-stop [PATH] | checkpoint [PATH] | threads N [blas M])"; return false; }
    std::string extra;
    if (is>>extra) { why="unexpected '"+extra+"' after the command"; return false; }
    return true;
}

FrontEndRunControl::FrontEndRunControl(std::string who, std::filesystem::path dir, std::string defaultSave)
    : itsWho(std::move(who)), itsDir(std::move(dir)), itsDefaultSave(std::move(defaultSave)), itsStarted(std::time(nullptr))
{
    namespace fs=std::filesystem;
    std::error_code ec;
    fs::create_directories(itsDir, ec);   // a control directory that cannot be made disables the file channel, never the run
    const std::string pid=std::to_string(::getpid());
    itsControl =itsDir/(pid+".control");
    itsProgress=itsDir/(pid+".progress");
    WriteHeartbeat("");
}

FrontEndRunControl::~FrontEndRunControl()
{
    std::error_code ec;
    std::filesystem::remove(itsProgress, ec);
}

void FrontEndRunControl::InstallSignalHandlers()
{
    struct sigaction once{};
    once.sa_handler=OnStopSignal;
    sigemptyset(&once.sa_mask);
    once.sa_flags=SA_RESETHAND;   // the SECOND SIGTERM / Ctrl-C takes the default action: it kills
    sigaction(SIGTERM, &once, nullptr);
    sigaction(SIGINT,  &once, nullptr);
    struct sigaction again{};
    again.sa_handler=OnStopSignal;
    sigemptyset(&again.sa_mask);
    sigaction(SIGUSR1, &again, nullptr);
}

void FrontEndRunControl::SetCurrent(const std::string& what)
{
    std::lock_guard<std::mutex> g(itsLock);
    itsCurrent=what;
}

std::vector<RunCommand> FrontEndRunControl::Poll()
{
    namespace fs=std::filesystem;
    std::lock_guard<std::mutex> g(itsLock);
    ++itsPolls;
    std::vector<RunCommand> cmds;
    std::string last;
    if (g_signalRequest.exchange(0)!=0)
    {
        RunCommand c; c.kind=RunCommand::Kind::SaveAndStop;
        cmds.push_back(c); last="save-and-stop (signal)";
    }
    std::error_code ec;
    if (fs::exists(itsControl, ec))
    {
        const fs::path taken=fs::path(itsControl.string()+".applying");
        fs::rename(itsControl, taken, ec);   // take the file first: a runctl write after this lands in a fresh one
        if (!ec)
        {
            std::ifstream in(taken);
            std::ofstream done(itsControl.string()+".done", std::ios::app);
            std::string line;
            while (std::getline(in, line))
            {
                if (Trim(line).empty() || Trim(line)[0]=='#') continue;
                RunCommand c; std::string why;
                if (ParseControlLine(line, c, why))
                {
                    cmds.push_back(c); last=Trim(line);
                    done<<Trim(line)<<"   # applied at poll "<<itsPolls<<"\n";
                }
                else
                {
                    std::cerr<<"[run control] IGNORED control line '"<<Trim(line)<<"': "<<why<<std::endl;
                    done<<Trim(line)<<"   # IGNORED: "<<why<<"\n";
                }
            }
            in.close();
            fs::remove(taken, ec);
        }
    }
    for (const RunCommand& c : cmds)
        if (c.kind==RunCommand::Kind::Stop || c.kind==RunCommand::Kind::SaveAndStop) itsStopped=true;
    if (itsStopped && cmds.empty()) { RunCommand c; c.kind=RunCommand::Kind::Stop; cmds.push_back(c); }   // sticky
    WriteHeartbeat(last);
    return cmds;
}

void FrontEndRunControl::WriteHeartbeat(const std::string& lastCommand)
{
    namespace fs=std::filesystem;
    const ThreadCounts t=CurrentThreads();
    std::ostringstream os;
    os<<"{\"pid\": "<<::getpid()<<", \"who\": "<<JsonString(itsWho)<<", \"current\": "<<JsonString(itsCurrent)
      <<", \"started\": "<<itsStarted<<", \"updated\": "<<std::time(nullptr)<<", \"polls\": "<<itsPolls
      <<", \"openmp\": "<<t.openmp<<", \"blas\": "<<t.blas<<", \"stopped\": "<<(itsStopped ? "true" : "false")
      <<", \"lastCommand\": "<<JsonString(lastCommand)<<", \"control\": "<<JsonString(itsControl.string())<<"}\n";
    std::error_code ec;
    const fs::path tmp=fs::path(itsProgress.string()+".tmp");
    { std::ofstream o(tmp, std::ios::trunc); if (!o) return; o<<os.str(); }
    fs::rename(tmp, itsProgress, ec);
}

} // namespace qchem::frontend
