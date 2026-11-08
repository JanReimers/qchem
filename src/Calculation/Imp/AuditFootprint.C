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

namespace
{
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
                                return {std::max(0.0, *limit-*used), "cgroup "+std::string(file)+" - memory.current ("+path+")"};
            if (path=="/") break;
            const auto slash=path.find_last_of('/');
            path = slash==0 ? "/" : path.substr(0, slash);
        }
    }
    if (const auto mi=ReadFile("/proc/meminfo"))
        if (const auto a=ParseMemAvailable(*mi)) return {*a, "MemAvailable (/proc/meminfo)"};
    throw std::runtime_error("QueryAvailableMemory: neither a cgroup v2 limit nor /proc/meminfo's MemAvailable could be read");
}

FootprintVerdict Judge(const AuditFootprint& f, double availableBytes)
{
    if (f.Low ()>availableBytes) return FootprintVerdict::Refuse;   // even the smallest it was ever measured at will not fit
    if (f.High()>availableBytes) return FootprintVerdict::Warn;     // it may not
    return FootprintVerdict::Fits;
}

report::json Describe(const AuditFootprint& f, const AvailableMemory& a)
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
        {"availableBytes", a.bytes}, {"availableSource", a.source}};
}

std::string ResourceRefusal::Summary() const
{
    auto GB=[](double b) {std::ostringstream os; os<<std::fixed<<std::setprecision(2)<<b/1.0e9<<" GB"; return os.str();};
    return "audit reference '"+footprint.reference+"' ("+when+") needs "+GB(footprint.Low())+" to "+GB(footprint.High())+" (estimate "
           +GB(footprint.estimate)+"), and "+GB(available.bytes)+" is available ["+available.source+"] -- "+escape;
}

} // namespace qchem
