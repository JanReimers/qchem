// File: Internal/LinearRadialMesh.C  Uniform radial mesh builder with trapezoidal r^2 weights.
// Elementary quadrature (NOT a transplanted special scheme): r_i evenly spaced on [start,stop],
// w_i = (trapezoid weight in r) * r_i^2, so sum_i w_i f(r_i) ~ integral_start^stop r^2 f dr.
module;
#include <cassert>
#include <utility>
module qchem.Mesh1.Radial;

namespace qcMesh1
{

RadialMesh LinearRadial(double start, double stop, int NumPoints)
{
    assert(NumPoints>=2);
    rvec_t R(NumPoints), W(NumPoints);
    double dr=(stop-start)/(NumPoints-1);
    for (int i=0; i<NumPoints; i++)
    {
        double r=start+i*dr;
        double trap=(i==0 || i==NumPoints-1) ? 0.5*dr : dr;
        R[i]=r;
        W[i]=trap*r*r;
    }
    return RadialMesh(std::move(R), std::move(W));
}

} //namespace qcMesh1
