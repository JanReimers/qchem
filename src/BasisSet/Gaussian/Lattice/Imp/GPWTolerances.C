// File: BasisSet/Gaussian/Lattice/Imp/GPWTolerances.C  See the interface.
module;
#include <sstream>
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Lattice.GPWTolerances;

namespace qchem::BasisSet::Gaussian
{

std::string GPWTolerances::Describe() const
{
    const GPWTolerances d;
    std::ostringstream os;
    if (vlocEps          != d.vlocEps)          os << " vlocEps="          << vlocEps;
    if (localPPRelCutoff != d.localPPRelCutoff) os << " localPPRelCutoff=" << localPPRelCutoff << " Ha";
    if (relFieldSharp    != d.relFieldSharp)    os << " relFieldSharp="    << relFieldSharp;
    if (screenEps        != d.screenEps)        os << " screenEps="        << screenEps;
    if (fieldSharp       != d.fieldSharp)       os << " fieldSharp="       << fieldSharp;
    if (relCutoff        != d.relCutoff)        os << " relCutoff="        << relCutoff << " Ha";
    if (densityEps       != d.densityEps)       os << " densityEps="       << densityEps;
    if (!mgridEcuts.empty()) { os << " mgridEcuts="; for (size_t i=0;i<mgridEcuts.size();++i) os << (i?",":"") << mgridEcuts[i]; }
    return os.str();
}

} // namespace qchem::BasisSet::Gaussian
