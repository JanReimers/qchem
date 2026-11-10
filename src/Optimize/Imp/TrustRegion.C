// File: Optimize/Imp/TrustRegion.C  The shared trust radius and its noise rules (the header says why).
module;
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
module qchem.Optimize.TrustRegion;

namespace qchem::Optimize
{

namespace
{
std::string Sci(double x) {char b[32]; std::snprintf(b, sizeof b, "%.3g", x); return b;}
}

TrustRegion::TrustRegion(const TrustOptions& options) : itsOpts(options), itsRadius(options.initial)
{
    if (!(itsOpts.initial>0.0 && itsOpts.min>0.0 && itsOpts.max>=itsOpts.initial))
        throw std::invalid_argument("TrustRegion: the radii must satisfy 0 < min, 0 < initial <= max");
}

Verdict TrustRegion::Judge(double change, double predicted, double noise, bool limited)
{
    if (change>noise)   // RULE 1
    {
        itsRadius*=0.25;
        if (itsRadius<itsOpts.min)
            return {Verdict::Kind::Collapsed, "the trust radius fell below "+Sci(itsOpts.min)+" after a value rise of "+Sci(change)
                                              +" (noise bound "+Sci(noise)+")"};
        return {Verdict::Kind::Reject, {}};
    }
    if (std::fabs(predicted)>noise)   // RULE 2
    {
        const double ratio = predicted<0.0 ? change/predicted : 0.0;
        if      (ratio>0.75 && limited) itsRadius=std::min(2.0*itsRadius, itsOpts.max);
        else if (ratio<0.25)            itsRadius=std::max(0.5*itsRadius, itsOpts.min);
    }
    return {Verdict::Kind::Accept, {}};
}

} // namespace qchem::Optimize
