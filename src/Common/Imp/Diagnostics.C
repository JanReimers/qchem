// File: Common/Imp/Diagnostics.C  The diagnostics registry (D-ENV step 2).  See the interface for the rules.
module;
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Diagnostics;

namespace qchem::Diagnostics
{

const std::vector<Entry>& Registry()
{
    // id, legacy env alias, what it prints.  NONE of these may change a computed number.
    static const std::vector<Entry> r={
        {"becke_atoms",      "GPW_BECKE_ATOMS",       "per-atom Becke share / kept / dropped breakdown"},
        {"becke_count",      "GPW_BECKE_COUNT",       "census of Becke per-point partition costs"},
        {"dm_rank",          "GPW_DM_RANK",           "rank / PSD census of D and the IPR of its Cholesky orbitals (the pin-21 canary)"},
        {"eh_trace",         "GPW_EH_TRACE",          "E_H G-space pairing trace"},
        {"field_spectrum",   "GPW_FIELD_SPECTRUM",    "cumulative |c|^2 spectrum of a fitted field"},
        {"gdm_trace",        "GPW_GDMTRACE",          "GDM line-search energy trace"},
        {"integrate_census", "GPW_INTEGRATE_CENSUS",  "NEW / REPEAT census of integrate-back calls"},
        {"kerker_spectrum",  "GPW_KERKER_SPECTRUM",   "Kerker residual power binned by G/G0"},
        {"localpp_timing",   "",                      "times the local-PP build"},   // no legacy alias: GPW_LOCALPP_RELCUTOFF is the numeric kappa knob
        {"mesh_ortho",       "GPW_MESH_ORTHO",        "plane-wave orthogonality error of the XC mesh, binned in |dG|; value = half-width of the dG box (default 1; 4 for a real onset)"},
        {"metal_trace",      "GPW_METALTRACE",        "shared-mu fill, per-block trace"},
        {"nl_per_l",         "GPW_NL_PER_L",          "bank the per-l non-local projector blocks (I0 diagnostic)"},
        {"phi_sparsity",     "GPW_PHI_SPARSITY",      "block-sparsity ceiling of the Phi table"},
        {"rho_negative",     "GPW_RHO_NEGATIVE",      "negative-rho census per XC route"},
        {"rss_trace",        "GPW_RSS_TRACE",         "resident-memory breadcrumbs in the lattice sum"},
        {"xc_alpha",         "GPW_XC_ALPHA",          "the XC mix's effective alpha each step"},
        {"xc_route",         "GPW_XCROUTE",           "which V_xc route fires each iteration"},
        {"angmesh_debug",    "QCHEM_ANGMESH_DEBUG",   "NNLS site-adapted angular-mesh debug print"},
        {"dump_h",           "QCHEM_DUMP_H",          "||F||, trace and max imaginary part of each Hamiltonian"},
        {"mom_scores",       "QCHEM_MOM_SCORES",      "sorted head of the MOM scores at each fill"},
        {"site_moments",     "QCHEM_SITE_MOMENTS",    "integrated site moments each iteration (to be promoted into the run report)"},
        {"u_trace",          "QCHEM_U_TRACE",         "per-refresh +U occupation line on stdout"},
    };
    return r;
}

static bool Known(const std::string& id)
{ for (const Entry& e : Registry()) if (e.id==id) return true; return false; }

std::map<std::string,std::string> ParseList(const std::string& list, std::vector<std::string>& unknown)
{
    std::map<std::string,std::string> out;
    std::stringstream ss(list);
    std::string tok;
    while (std::getline(ss, tok, ','))
    {
        if (tok.empty()) continue;
        const size_t eq=tok.find('=');
        const std::string id=tok.substr(0, eq), val=(eq==std::string::npos) ? "1" : tok.substr(eq+1);
        if (id=="list") { out["list"]="1"; continue; }
        if (!Known(id)) { unknown.push_back(id); continue; }
        out[id]=val;
    }
    return out;
}

void Describe(std::ostream& os)
{
    os<<"[diagnostics] QCHEM_DIAGNOSTICS=<id>[=value],...   (an old GPW_*/QCHEM_* name still works as an alias):\n";
    for (const Entry& e : Registry()) os<<"  "<<e.id<<std::string(std::max<size_t>(1, 18-e.id.size()),' ')<<e.legacyEnv
                                          <<std::string(std::max<size_t>(1, 24-e.legacyEnv.size()),' ')<<e.what<<"\n";
}

namespace
{
std::mutex& Mu() { static std::mutex m; return m; }
std::map<std::string,std::string>& Overrides() { static std::map<std::string,std::string> o; return o; }   // Scoped (tests)
std::map<std::string,std::string>& FromEnv()
{
    static std::map<std::string,std::string> m=[]
    {
        std::vector<std::string> unknown;
        std::map<std::string,std::string> parsed;
        if (const char* s=std::getenv("QCHEM_DIAGNOSTICS")) parsed=ParseList(s, unknown);
        for (const std::string& u : unknown)
        {
            std::string near;
            for (const Entry& e : Registry())                       // the typo catch: ids sharing a prefix or substring
                if (e.id.find(u)!=std::string::npos || u.find(e.id)!=std::string::npos
                    || (u.size()>=3 && e.id.compare(0,3,u,0,3)==0)) near+=" "+e.id;
            std::cerr<<"[diagnostics] WARNING: unknown id '"<<u<<"' in QCHEM_DIAGNOSTICS"
                     <<(near.empty() ? std::string(" (QCHEM_DIAGNOSTICS=list shows them all)") : "; did you mean:"+near)<<std::endl;
        }
        if (parsed.count("list")) { Describe(std::cout); parsed.erase("list"); }
        for (const Entry& e : Registry())                           // legacy alias: ON unless "0"; its value is the argument
            if (e.legacyEnv.empty()) continue;
            else if (const char* v=std::getenv(e.legacyEnv.c_str()); v && std::string(v)!="0" && !parsed.count(e.id))
                parsed[e.id]=*v ? v : "1";
        return parsed;
    }();
    return m;
}
}

std::optional<std::string> Value(const std::string& id)
{
    if (!Known(id)) throw std::logic_error("Diagnostics: '"+id+"' is not a registered diagnostic (add it to Registry())");
    std::lock_guard<std::mutex> lk(Mu());
    if (auto it=Overrides().find(id); it!=Overrides().end()) return it->second=="off" ? std::nullopt : std::optional<std::string>(it->second);
    if (auto it=FromEnv().find(id); it!=FromEnv().end()) return it->second;
    return std::nullopt;
}

bool Enabled(const std::string& id) { return Value(id).has_value(); }

Scoped::Scoped(const std::string& id, const std::string& value) : itsId(id)
{
    if (!Known(id)) throw std::logic_error("Diagnostics::Scoped: '"+id+"' is not registered");
    std::lock_guard<std::mutex> lk(Mu());
    if (auto it=Overrides().find(id); it!=Overrides().end()) itsPrev=it->second;
    Overrides()[id]=value;
}
Scoped::~Scoped()
{
    std::lock_guard<std::mutex> lk(Mu());
    if (itsPrev) Overrides()[itsId]=*itsPrev; else Overrides().erase(itsId);
}

} // namespace qchem::Diagnostics
