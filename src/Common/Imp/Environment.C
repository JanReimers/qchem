// File: Common/Imp/Environment.C  See the interface.
module;
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <set>
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
}
