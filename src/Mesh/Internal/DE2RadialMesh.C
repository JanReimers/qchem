// File: Internal/DE2RadialMesh.C  The DE2 double-exponential radial mesh builder (Mitani 2011, Mitani & Yoshioka 2012;
// Dasgupta & Herbert 2017 eqs. 6-9, where it is the radial rule of SG-2 / SG-3).
//
// The map r = t(x) = exp(alpha x - e^{-x}) takes x in (-inf, inf) onto r in (0, inf).  Under it the radial integral becomes
//     integral_0^inf F(r) r^2 dr = integral_-inf^inf F(t(x)) t(x)^2 t'(x) dx,   t'(x) = t(x) (alpha + e^{-x}),
// and the trapezoid rule on equally spaced x_i = x_0 + i h gives the nodes and weights (r^2 jacobian folded in, as for MHL)
//     r_i = exp(alpha x_i - e^{-x_i}),    w_i = h exp(3 alpha x_i - 3 e^{-x_i}) (alpha + e^{-x_i}) = h r_i^3 (alpha + e^{-x_i}).
// The transformed integrand decays DOUBLE exponentially at both ends: as x -> -inf, r -> 0 like exp(-e^{-x}); as x -> +inf, r grows
// like e^{alpha x} and an atom's density decays like exp(-c r) or exp(-a r^2) of that.  For an integrand analytic in a strip the
// trapezoid error then falls like exp(-c N / log N) (Takahasi-Mori; Mitani 2011 section 3), so a smooth pseudo-density gains most.
//
// THE TRUNCATION is Mitani's: the N nodes span [rMin, rMax] (his rMin = 1e-7, rMax = 10 R with R the outer valence orbital radius),
// both ends included, so h = (x(rMax) - x(rMin)) / (N - 1).  The weights at both ends are the double-exponentially small tails
// the infinite sum drops.  alpha is a SHAPE parameter (how the nodes are shared between core and tail): it is a calibration, stated
// by the caller and never folded into this formula (qcMesh::kDE2AlphaCalibrated in qchem.Mesh holds the measured value and its record).
module;
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
module qchem.Mesh.Radial;

namespace qchem::qcMesh
{

//! x with t(x) = r: alpha x - e^{-x} = ln r.  The left side increases strictly (its slope alpha + e^{-x} > 0), so bisection on a
//! bracket that holds every double r > 0 converges; 200 halvings reach the bracket's roundoff.
static double DE2Abscissa(double r, double alpha)
{
    const double lnr=std::log(r);
    double lo=-40.0, hi=800.0/alpha;   // t(-40) = exp(-e^40) underflows to 0; t(800/alpha) = e^800 overflows
    for (int it=0; it<200 && hi-lo>1e-15*std::max(1.0,std::fabs(lo)); it++)
    {
        const double mid=0.5*(lo+hi);
        (alpha*mid-std::exp(-mid) < lnr ? lo : hi)=mid;
    }
    return 0.5*(lo+hi);
}

RadialMesh DE2Radial(int NumPoints, double alpha, double rMin, double rMax)
{
    if (NumPoints<2) throw std::invalid_argument("DE2Radial: needs at least 2 points, got "+std::to_string(NumPoints));
    if (!(alpha>0.0)) throw std::invalid_argument("DE2Radial: the mapping parameter alpha must be positive, got "+std::to_string(alpha));
    if (!(rMin>0.0 && rMax>rMin))
        throw std::invalid_argument("DE2Radial: needs 0 < rMin < rMax, got rMin="+std::to_string(rMin)+" rMax="+std::to_string(rMax));
    const double x0=DE2Abscissa(rMin,alpha), x1=DE2Abscissa(rMax,alpha), h=(x1-x0)/(NumPoints-1);
    rvec_t R(NumPoints), W(NumPoints);
    for (int i=0; i<NumPoints; i++)
    {
        const double x=x0+i*h, e=std::exp(-x), r=std::exp(alpha*x-e);
        R[i]=r;
        W[i]=h*r*r*r*(alpha+e);
    }
    return RadialMesh(std::move(R), std::move(W));
}

} //namespace qchem::qcMesh
