// File: Calculation/Imp/RasterPilotXC.C  The uniform XC term's L0 demand (the interface says why).
module;
#include <algorithm>
#include <optional>
module qchem.RasterPilotXC;
import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // GridAccuracy::MeasuredDemand / MeasuredPowerDemand

namespace qchem
{

UniformXCDemand::Verdict UniformXCDemand::Round(double xL0, double d, double target)
{
    using GA=BasisSet::Gaussian::GridAccuracy;
    Verdict v;
    const bool moved=itsBefore && xL0>itsBefore->xL0*(1.0+1e-9);
    if (moved && !(d<kImprovedRatio*itsBefore->d) && !itsRefuted) { itsRefuted=true; v.refuted=true; }
    v.over=d>target;
    if (!itsRefuted)
    {
        double need=GA::MeasuredDemand(d, xL0, target);   // 0 = meets the target here
        if (need>0.0 && moved)
            if (const double p=GA::MeasuredPowerDemand(itsBefore->xL0, itsBefore->d, xL0, d, target); p>0.0) need=p;
        v.needL0 = need>0.0 ? std::min(need, kMaxStep*xL0) : 0.0;
    }
    itsStalled=v.over && itsRefuted;
    itsBefore=Point{xL0, d};
    return v;
}

} // namespace qchem
