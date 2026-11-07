// File: BasisSet/Gaussian/Lattice/Imp/GPWTolerances.C  See the interface.
module;
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
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
    if (VlocEpsDerived())                       os << " (derived: rigorous " << vlocEpsRigorous << " / eta_routing " << vlocRoutingEta << ")";
    if (localPPRelCutoff != d.localPPRelCutoff) os << " localPPRelCutoff=" << localPPRelCutoff << " Ha";
    if (relFieldSharp    != d.relFieldSharp)    os << " relFieldSharp="    << relFieldSharp;
    if (screenEps        != d.screenEps)        os << " screenEps="        << screenEps;
    if (fieldSharp       != d.fieldSharp)       os << " fieldSharp="       << fieldSharp;
    if (relCutoff        != d.relCutoff)        os << " relCutoff="        << relCutoff << " Ha";
    if (densityEps       != d.densityEps)       os << " densityEps="       << densityEps;
    if (rungGateC        != d.rungGateC)        os << " rungGateC="        << rungGateC;
    if (relCutoffSafety  != d.relCutoffSafety)  os << " relCutoffSafety="  << relCutoffSafety;
    if (rasterX          != d.rasterX)          os << " rasterX="          << rasterX;
    if (vlocRasterX      != d.vlocRasterX)      os << " vlocRasterX="      << vlocRasterX;
    if (xcRasterX        != d.xcRasterX)        os << " xcRasterX="        << xcRasterX;
    if (!mgridEcuts.empty()) { os << " mgridEcuts="; for (size_t i=0;i<mgridEcuts.size();++i) os << (i?",":"") << mgridEcuts[i]; }
    return os.str();
}

double AliasTail::Bound(double x) const {return S*std::pow(x,n)*std::exp(-x);}

double AliasTail::Exponent(double eps) const
{
    // f(L) = ln S + n ln L - L is decreasing for L > n; its maximum is at L = n.  Newton on f(L) = ln eps from the bare guess.
    if (!(eps>0.0)) throw std::invalid_argument("AliasTail::Exponent: eps must be > 0");
    const double lnEps=std::log(eps);
    if (std::log(S)+n*std::log(n)-n <= lnEps) throw std::invalid_argument("AliasTail::Exponent: eps exceeds the tail's maximum -- no raster needed");
    double L=std::max(n+1.0, std::log(S/eps));
    for (int it=0; it<100; ++it)
    {
        const double f=std::log(S)+n*std::log(L)-L-lnEps, df=n/L-1.0;
        const double dL=-f/df;
        L+=dL;
        if (L<=n) L=n*(1.0+1e-12)+1e-12;
        if (std::fabs(dL)<1e-13*L) return L;
    }
    throw std::logic_error("AliasTail::Exponent: Newton did not converge");
}

double GridAccuracy::Calibrated(const AliasTail& t, double etaTerm) const
{
    if (!(etaTerm>0.0 && etaTerm<=1.0)) throw std::invalid_argument("GridAccuracy: eta = delivered/bound must be in (0, 1] (1 = rigorous)");
    const double target=eps/etaTerm;
    return t.Bound(t.n) <= target ? t.n : t.Exponent(target);
}

GPWTolerances GridAccuracy::Stated(GPWTolerances t) const
{
    if (!eta) throw std::logic_error("GridAccuracy::Stated: no stated calibration -- this run self-calibrates (the pilot audit)");
    t.rasterX    =Calibrated(kDensityTail, eta->density);
    t.vlocRasterX=Calibrated(kVlocLongTail, eta->vlocLong);
    t.xcRasterX  =t.rasterX;
    return t;
}

double GridAccuracy::MeasuredDemand(double delivered, double xRun, double target)
{
    if (!(target>0.0)) throw std::invalid_argument("GridAccuracy::MeasuredDemand: the target must be > 0");
    if (!(delivered>target)) return 0.0;
    if (!(xRun>0.0)) throw std::invalid_argument("GridAccuracy::MeasuredDemand: a term over its target on a raster with no alias exponent");
    const AliasTail measured{delivered/(std::pow(xRun,kAliasTailPower)*std::exp(-xRun)), kAliasTailPower};
    return measured.Bound(measured.n) <= target ? measured.n : measured.Exponent(target);
}

double GridAccuracy::MeasuredPowerDemand(double x1, double d1, double x2, double d2, double target)
{
    if (!(target>0.0)) throw std::invalid_argument("GridAccuracy::MeasuredPowerDemand: the target must be > 0");
    if (!(d2>target)) return 0.0;
    if (!(x2>x1 && x1>0.0 && d1>d2 && d2>0.0)) return 0.0;   // no decrease measured: no slope to extrapolate (the caller falls back)
    const double k=std::log(d1/d2)/std::log(x2/x1);         // d ~ x^{-k}
    return x2*std::pow(d2/target, 1.0/k);
}

std::string GridAccuracy::Describe(const GPWTolerances& used) const
{
    std::ostringstream os;
    if (!(eps>0.0))
    {
        os<<"NO target (eps 0): the rasters as stated, x(density, V_long, XC L0)=("<<used.rasterX<<", "<<used.vlocRasterX<<", "<<used.xcRasterX
          <<") (0 = the ball's own raster)";
        return os.str();
    }
    os<<std::setprecision(3)<<"eps="<<eps<<" Ha";
    if (eta) os<<"  eta(V_long, density)=("<<eta->vlocLong<<", "<<eta->density<<") STATED";
    else     os<<"  eta SELF-CALIBRATED (the pilot audit sizes the rasters)";
    os<<"  rigorous x=G_N^2/p: V_long "<<std::setprecision(4)<<Rigorous(kVlocLongTail)<<", density "<<Rigorous(kDensityTail)
      <<";  XC raster: MEASURED (no rigorous bound)"
      <<";  USED: raster x density "<<used.rasterX<<", V_long "<<used.vlocRasterX<<", XC L0 "<<used.xcRasterX<<" (0 = the ball's own raster)"
      <<", vlocEps="<<std::setprecision(3)<<used.vlocEps;
    if (used.VlocEpsDerived()) os<<" (derived: rigorous per-pair "<<used.vlocEpsRigorous<<" / calibration eta_routing "<<used.vlocRoutingEta<<")";
    else                     os<<(used.vlocEps>0.0 ? " (STATED)" : " (no-target fallback)");
    os
      <<", kappa="<<used.localPPRelCutoff<<" Ha (L "<<used.localPPRelCutoff/2<<"), screenEps="<<used.screenEps<<", densityEps="<<used.densityEps;
    return os.str();
}

} // namespace qchem::BasisSet::Gaussian
