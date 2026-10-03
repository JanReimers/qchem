// File: BasisSet/Gaussian/Lattice/Imp/GPWTolerances.C  See the interface.
module;
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Lattice.GPWTolerances;
import qchem.Environment;

namespace qchem::BasisSet::Gaussian
{

std::string GPWTolerances::Describe() const
{
    const GPWTolerances d;
    std::ostringstream os;
    if (vlocEps          != d.vlocEps)          os << " vlocEps="          << vlocEps;
    if (localPPRelCutoff != d.localPPRelCutoff) os << " localPPRelCutoff=" << localPPRelCutoff << " Ha";
    if (relFieldSharp    != d.relFieldSharp)    os << " relFieldSharp="    << relFieldSharp;
    if (!mgridEcuts.empty()) { os << " mgridEcuts="; for (size_t i=0;i<mgridEcuts.size();++i) os << (i?",":"") << mgridEcuts[i]; }
    return os.str();
}

std::string ApplyEnvOverrides(GPWTolerances& t)
{
    std::ostringstream said;
    auto num=[&](const char* name, double& v)
    {
        if (const char* s=std::getenv(name)) { v=std::atof(s); said << " " << name << "=" << s; }
    };
    num("GPW_VLOC_EPS",          t.vlocEps);
    num("GPW_LOCALPP_RELCUTOFF", t.localPPRelCutoff);
    num("GPW_RELFIELDSHARP",     t.relFieldSharp);
    if (const char* s=std::getenv("GPW_MGRID_ECUTS"))
    {
        t.mgridEcuts.clear();
        for (std::string list(s); !list.empty();)
        {
            const size_t c=list.find(',');
            t.mgridEcuts.push_back(std::atof(list.substr(0,c).c_str()));
            list = c==std::string::npos ? std::string() : list.substr(c+1);
        }
        said << " GPW_MGRID_ECUTS=" << s;
    }
    return said.str();
}

} // namespace qchem::BasisSet::Gaussian
