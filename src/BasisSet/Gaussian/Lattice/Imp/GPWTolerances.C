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
    if (localPPRelCutoff != d.localPPRelCutoff) os << " localPPRelCutoff=" << localPPRelCutoff << " Ha";
    if (relFieldSharp    != d.relFieldSharp)    os << " relFieldSharp="    << relFieldSharp;
    if (screenEps        != d.screenEps)        os << " screenEps="        << screenEps;
    if (fieldSharp       != d.fieldSharp)       os << " fieldSharp="       << fieldSharp;
    if (relCutoff        != d.relCutoff)        os << " relCutoff="        << relCutoff << " Ha";
    if (densityEps       != d.densityEps)       os << " densityEps="       << densityEps;
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

std::string GridAccuracy::Describe(const GPWTolerances& used) const
{
    std::ostringstream os;
    os<<std::setprecision(3)<<"eps="<<eps<<" Ha  eta(V_long, density)=("<<eta.vlocLong<<", "<<eta.density<<")"
      <<"  raster x=G_N^2/p demanded (rigorous / calibrated): V_long "<<std::setprecision(4)<<Rigorous(kVlocLongTail)<<" / "<<Calibrated(kVlocLongTail,eta.vlocLong)
      <<", density "<<Rigorous(kDensityTail)<<" / "<<Calibrated(kDensityTail,eta.density)
      <<";  XC raster: MEASURED (no rigorous bound)"
      <<";  USED today (step 3 derives them): raster = the RasterPolicy, vlocEps="<<std::setprecision(3)<<used.vlocEps
      <<", kappa="<<used.localPPRelCutoff<<" Ha (L "<<used.localPPRelCutoff/2<<"), screenEps="<<used.screenEps<<", densityEps="<<used.densityEps;
    return os.str();
}

} // namespace qchem::BasisSet::Gaussian
