// File: Calculation/Imp/RasterPilotXC.C  The uniform XC term's two-arm demand (the interface says why).
module;
#include <algorithm>
#include <optional>
module qchem.RasterPilotXC;
import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // GridAccuracy::MeasuredDemand / MeasuredPowerDemand

namespace qchem
{

UniformXCDemand::Verdict UniformXCDemand::Round(double xRho, double xL0, double d, double target, bool belowRho)
{
    using GA=BasisSet::Gaussian::GridAccuracy;
    Verdict v;
    bool movedRho=false, movedL0=false;
    if (itsBefore)
    {
        movedRho=xRho>itsBefore->xRho*(1.0+1e-9);
        movedL0 =xL0 >itsBefore->xL0 *(1.0+1e-9);
        if (!(d<kImprovedRatio*itsBefore->d))   // the error did not fall: an arm that moved ALONE did not cause it
        {
            if (movedRho && !movedL0 && !itsRhoRefuted) { itsRhoRefuted=true; v.refuted="density raster"; }
            if (movedL0 && !movedRho && !itsL0Refuted)  { itsL0Refuted=true;  v.refuted="L0 raster"; }
        }
    }
    // One arm's demand on its own axis: the power law once that axis has moved since the last round, else the rigorous shape's MeasuredDemand.
    auto armDemand=[&](double xNow, std::optional<double> xThen, bool moved)
    {
        double need=GA::MeasuredDemand(d, xNow, target);   // 0 = meets the target here
        if (need>0.0 && xThen && moved)
            if (const double p=GA::MeasuredPowerDemand(*xThen, itsBefore->d, xNow, d, target); p>0.0) need=p;
        return need>0.0 ? std::min(need, kMaxStep*xNow) : 0.0;
    };
    v.rhoArm=belowRho && !itsRhoRefuted;
    v.needRho=v.rhoArm        ? armDemand(xRho, itsBefore ? std::optional<double>(itsBefore->xRho) : std::nullopt, movedRho) : 0.0;
    v.needL0 =!itsL0Refuted   ? armDemand(xL0,  itsBefore ? std::optional<double>(itsBefore->xL0)  : std::nullopt, movedL0)  : 0.0;
    v.over=d>target;
    v.unanswered=v.over && v.needRho==0.0 && v.needL0==0.0;
    itsBefore=Point{xRho, xL0, d};
    return v;
}

} // namespace qchem
