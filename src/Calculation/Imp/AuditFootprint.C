// File: Calculation/Imp/AuditFootprint.C  The audit reference's byte count, its calibrated estimate, and the memory the run has.
module;
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
module qchem.AuditFootprint;

namespace qchem
{

AuditBytes CountAuditBytes(const AuditReferenceSize& s)
{
    double sumS=0.0, sumNS=0.0;   // Sum_b s_b and Sum_b nf_b s_b
    for (const AuditBlock& b : s.blocks)
    {
        const double bytes = b.real ? 8.0 : 16.0;
        sumS +=bytes;
        sumNS+=bytes*double(b.nFunctions);
    }
    AuditBytes a;
    a.phiTables   =s.meshPoints*sumNS;
    a.rasterFields=s.rasterPoints*sumS*double(s.spinChannels);
    return a;
}

AuditFootprint EstimateAuditFootprint(const AuditReferenceSize& s, const AuditFootprintCalibration& cal)
{
    AuditFootprint f;
    f.reference=s.name;
    f.size     =s;
    f.bytes    =CountAuditBytes(s);
    f.cal      =cal;
    f.estimate =cal.eta*f.bytes.Known() + cal.kappa*s.rasterPoints;   // THE ONE PLACE the byte count and the calibration combine
    return f;
}

std::optional<double> ParseCgroupLimit(const std::string& text)
{
    std::istringstream is(text);
    std::string w; is>>w;
    if (w=="max") return std::nullopt;
    std::size_t used=0;
    double v=0.0;
    try { v=std::stod(w, &used); } catch (const std::exception&) { used=0; }
    if (w.empty() || used!=w.size()) throw std::runtime_error("ParseCgroupLimit: '"+text+"' is neither a byte count nor 'max'");
    return v;
}

std::optional<double> ParseMemAvailable(const std::string& meminfo)
{
    std::istringstream is(meminfo);
    std::string line;
    while (std::getline(is, line))
        if (line.rfind("MemAvailable:", 0)==0)
        {
            std::istringstream ls(line.substr(13));
            double kB=0.0; ls>>kB;
            if (!ls) throw std::runtime_error("ParseMemAvailable: unreadable line '"+line+"'");
            return kB*1024.0;
        }
    return std::nullopt;
}

std::optional<double> ParseCgroupReclaimable(const std::string& memoryStat)
{
    std::istringstream is(memoryStat);
    std::string key;
    double v=0.0, sum=0.0;
    bool found=false;
    while (is>>key>>v)
        if (key=="inactive_file" || key=="slab_reclaimable") { sum+=v; found=true; }
    if (found) return sum;
    return std::nullopt;
}

std::optional<double> ParseOwnResidentNonFile(const std::string& procStatus)
{
    std::istringstream is(procStatus);
    std::string line;
    double kB=0.0;
    bool found=false;
    while (std::getline(is, line))
        for (const char* key : {"RssAnon:", "RssShmem:"})
            if (line.rfind(key, 0)==0)
            {
                std::istringstream ls(line.substr(std::string(key).size()));
                double v=0.0; ls>>v;
                if (!ls) throw std::runtime_error("ParseOwnResidentNonFile: unreadable line '"+line+"'");
                kB+=v; found=true;
            }
    if (found) return kB*1024.0;
    return std::nullopt;
}

namespace
{
std::string GB(double b) {std::ostringstream os; os<<std::fixed<<std::setprecision(2)<<b/1.0e9<<" GB"; return os.str();}

std::optional<std::string> ReadFile(const std::string& path)
{
    std::ifstream f(path);
    if (!f) return std::nullopt;
    std::ostringstream os; os<<f.rdbuf();
    return os.str();
}
} // namespace

AvailableMemory QueryAvailableMemory()
{
    // cgroup v2: /proc/self/cgroup reads "0::/<path>"; the limits sit in /sys/fs/cgroup/<path>/memory.{high,max}.  The NEAREST limit up
    // the tree binds (memsafe sets MemoryHigh on its own scope), measured against that cgroup's own usage.
    if (const auto cg=ReadFile("/proc/self/cgroup"))
    {
        std::istringstream is(*cg);
        std::string line, path;
        while (std::getline(is, line)) if (line.rfind("0::", 0)==0) path=line.substr(3);
        while (!path.empty())
        {
            const std::string dir="/sys/fs/cgroup"+(path=="/" ? std::string() : path);
            for (const char* file : {"memory.high", "memory.max"})
                if (const auto t=ReadFile(dir+"/"+file))
                    if (const auto limit=ParseCgroupLimit(*t))
                        if (const auto cur=ReadFile(dir+"/memory.current"))
                            if (const auto used=ParseCgroupLimit(*cur))
                            {
                                // The page cache and reclaimable slab the scope holds are not commitments (the kernel drops them before it
                                // throttles): count only the rest as used.  Without a readable memory.stat the whole of memory.current is
                                // counted, and the source says so.
                                const auto st =ReadFile(dir+"/memory.stat");
                                const auto rec=st ? ParseCgroupReclaimable(*st) : std::nullopt;
                                const std::string head="cgroup "+std::string(file);
                                if (!rec)
                                    return {std::max(0.0, *limit-*used), head+" - memory.current ("+path+"; no memory.stat: page cache counted as used)"};
                                const double nonReclaimable=std::max(0.0, *used-*rec);
                                return {std::max(0.0, *limit-nonReclaimable),
                                        head+" - (memory.current - reclaimable) ("+path+"): current "+GB(*used)+", reclaimable "+GB(*rec)
                                        +" (inactive_file + slab_reclaimable)"};
                            }
            if (path=="/") break;
            const auto slash=path.find_last_of('/');
            path = slash==0 ? "/" : path.substr(0, slash);
        }
    }
    if (const auto mi=ReadFile("/proc/meminfo"))
        if (const auto a=ParseMemAvailable(*mi)) return {*a, "MemAvailable (/proc/meminfo)"};
    throw std::runtime_error("QueryAvailableMemory: neither a cgroup v2 limit nor /proc/meminfo's MemAvailable could be read");
}

std::optional<double> ParseMemTotal(const std::string& meminfo)
{
    std::istringstream is(meminfo);
    std::string line;
    while (std::getline(is, line))
        if (line.rfind("MemTotal:", 0)==0)
        {
            std::istringstream ls(line.substr(9));
            double kB=0.0; ls>>kB;
            if (!ls) throw std::runtime_error("ParseMemTotal: unreadable line '"+line+"'");
            return kB*1024.0;
        }
    return std::nullopt;
}

MemoryBudget QueryMemoryBudget()
{
    // The same walk as QueryAvailableMemory, but the LIMIT itself: the nearest memory.high / memory.max up the tree binds.
    if (const auto cg=ReadFile("/proc/self/cgroup"))
    {
        std::istringstream is(*cg);
        std::string line, path;
        while (std::getline(is, line)) if (line.rfind("0::", 0)==0) path=line.substr(3);
        while (!path.empty())
        {
            const std::string dir="/sys/fs/cgroup"+(path=="/" ? std::string() : path);
            for (const char* file : {"memory.high", "memory.max"})
                if (const auto t=ReadFile(dir+"/"+file))
                    if (const auto limit=ParseCgroupLimit(*t))
                        return {*limit, "cgroup "+std::string(file)+" ("+path+")"};
            if (path=="/") break;
            const auto slash=path.find_last_of('/');
            path = slash==0 ? "/" : path.substr(0, slash);
        }
    }
    if (const auto mi=ReadFile("/proc/meminfo"))
        if (const auto t=ParseMemTotal(*mi)) return {*t, "MemTotal (/proc/meminfo)"};
    throw std::runtime_error("QueryMemoryBudget: neither a cgroup v2 limit nor /proc/meminfo's MemTotal could be read");
}

ResourceRoom QueryResourceRoom()
{
    ResourceRoom r;
    r.limit=QueryMemoryBudget();
    r.live =QueryAvailableMemory();
    const auto st=ReadFile("/proc/self/status");
    const auto own=st ? ParseOwnResidentNonFile(*st) : std::nullopt;
    r.ownBytes =own ? *own : 0.0;
    r.ownSource=own ? "RssAnon + RssShmem (/proc/self/status)" : "/proc/self/status unreadable: 0 counted";
    return r;
}

std::string DescribeRoom(const ResourceRoom& r)
{
    return "fixed room "+GB(r.Fixed())+" ["+r.limit.source+" "+GB(r.limit.bytes)+" less this run's own "+GB(r.ownBytes)+", "+r.ownSource
           +"]; live "+GB(r.live.bytes)+" ["+r.live.source+"]";
}

FootprintVerdict Judge(const AuditFootprint& f, const ResourceRoom& r)
{
    if (f.Low ()>r.Fixed()) return FootprintVerdict::Refuse;   // even the smallest it was ever measured at will not fit in the fixed room
    if (f.High()>r.Warn ()) return FootprintVerdict::Warn;     // it may not, in the fixed room or beside what the scope holds now
    return FootprintVerdict::Fits;
}

report::json Describe(const AuditFootprint& f, const ResourceRoom& r)
{
    const AuditReferenceSize& s=f.size;
    std::size_t nReal=0; for (const auto& b : s.blocks) if (b.real) ++nReal;
    return report::json{
        {"reference", f.reference},
        {"sizes", {{"rasterPoints", s.rasterPoints}, {"meshPoints", s.meshPoints}, {"blocks", s.blocks.size()}, {"realBlocks", nReal},
                   {"spinChannels", s.spinChannels}}},
        {"byteCount", {{"phiTables", f.bytes.phiTables}, {"rasterFields", f.bytes.rasterFields}, {"known", f.bytes.Known()}}},
        {"calibration", {{"eta", f.cal.eta}, {"kappaBytesPerPoint", f.cal.kappa}, {"record", "kAuditFootprintRecord (20 audits, 4 runs)"}}},
        {"band", {{"low", f.cal.bandLow}, {"high", f.cal.bandHigh}}},
        {"estimateBytes", f.estimate}, {"lowBytes", f.Low()}, {"highBytes", f.High()},
        {"fixedRoomBytes", r.Fixed()}, {"limitBytes", r.limit.bytes}, {"limitSource", r.limit.source},
        {"ownBytes", r.ownBytes}, {"ownSource", r.ownSource},
        {"liveBytes", r.live.bytes}, {"liveSource", r.live.source}};
}

std::string ResourceRefusal::Summary() const
{
    return "audit reference '"+footprint.reference+"' ("+when+") needs "+GB(footprint.Low())+" to "+GB(footprint.High())+" (estimate "
           +GB(footprint.estimate)+"); "+DescribeRoom(room)+" -- "+escape;
}

} // namespace qchem
