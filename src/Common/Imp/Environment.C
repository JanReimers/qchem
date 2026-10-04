// File: Common/Imp/Environment.C  See the interface.
module;
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <set>
#include <vector>
#include <string>
module qchem.Environment;

namespace qchem
{
const char* Env(const char* name, const char* legacy)
{
    if (const char* v=std::getenv(name)) return v;
    if (!legacy) return nullptr;
    const char* v=std::getenv(legacy);
    if (v)
    {
        static std::mutex mu; static std::set<std::string> told;
        std::lock_guard<std::mutex> lk(mu);
        if (told.insert(legacy).second)
            std::cerr<<"[env] "<<legacy<<" is DEPRECATED (it is not specific to the Gaussian-plane-wave basis); use "<<name<<std::endl;
    }
    return v;
}

std::vector<RetiredVariable> RetiredEnvironmentSet()
{
    // name -> the deck key that replaces it.  D-ENV step 6a: the GPW tolerances and the Becke XC recipe.  Each Becke name also had a
    // GPW_BECKE_* deprecated alias, retired with it.
    static const std::vector<RetiredVariable> kRetired={
        {"GPW_VLOC_EPS","solid.tolerances.vlocEps"},        {"GPW_LOCALPP_RELCUTOFF","solid.tolerances.localPPRelCutoff"},
        {"GPW_RELFIELDSHARP","solid.tolerances.relFieldSharp"}, {"GPW_MGRID_ECUTS","solid.tolerances.mgridEcuts"},
        {"GPW_SCREEN_EPS","solid.tolerances.screenEps"},    {"GPW_FIELDSHARP","solid.tolerances.fieldSharp"},
        {"GPW_RELCUTOFF","solid.tolerances.relCutoff"},     {"GPW_DENSITY_EPS","solid.tolerances.densityEps"},
        {"QCHEM_BECKE_NR","solid.xcMesh.nRadial"},          {"GPW_BECKE_NR","solid.xcMesh.nRadial"},
        {"QCHEM_BECKE_ALPHA","solid.xcMesh.mhlAlpha"},      {"GPW_BECKE_ALPHA","solid.xcMesh.mhlAlpha"},
        {"QCHEM_BECKE_L","solid.xcMesh.angularDegree"},     {"GPW_BECKE_L","solid.xcMesh.angularDegree"},
        {"QCHEM_BECKE_ROT","solid.xcMesh.angRot"},          {"GPW_BECKE_ROT","solid.xcMesh.angRot"},
        {"QCHEM_BECKE_EPS","solid.xcMesh.beckeEps"},        {"GPW_BECKE_EPS","solid.xcMesh.beckeEps"},
        // 6b: the declared CP2K deviations (RunPolicy) are the deck's `solid.policy` block
        {"CP2K_COMPAT","solid.policy.cp2kCompat"},          {"QCHEM_DM_LOWRANK","solid.policy.dmLowRank"},
        {"GPW_STREAM_FOLD","solid.policy.streamFold"},      {"QCHEM_MIX_RHO_M","solid.policy.mixRhoM"},
        {"QCHEM_XC_DM_SOURCE","solid.policy.xcFromDM"},     {"GPW_XC_DM_SOURCE","solid.policy.xcFromDM"},
        {"QCHEM_IMPOSE_SYMMETRY","solid.policy.imposeSymmetry"}, {"QCHEM_BECKE_XC","solid.policy.beckeXC"},
        {"GPW_DAWARE_SCREEN","solid.policy.dAwareScreen"},  {"QCHEM_U_EIGEN","solid.policy.hubbardEigen"}};
    std::vector<RetiredVariable> set;
    for (const auto& r : kRetired) if (std::getenv(r.name.c_str())) set.push_back(r);
    return set;
}

void WarnRetiredEnvironment()
{
    static std::once_flag once;
    std::call_once(once, []
    {
        for (const auto& r : RetiredEnvironmentSet())
            std::cerr<<"[env] "<<r.name<<" is RETIRED and IGNORED: the input deck is the only way to set it (deck key '"<<r.deckKey
                     <<"', or --set "<<r.deckKey<<"=<value>)"<<std::endl;
    });
}
}
