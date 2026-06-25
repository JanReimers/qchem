// File: Internal/LinearRadialMesh.C  Uniform radial mesh with trapezoidal r^2 weights.
// Elementary quadrature (NOT a transplanted special scheme): r_i evenly spaced on [start,stop],
// w_i = (trapezoid weight in r) * r_i^2, so sum_i w_i f(r_i) ~ integral_start^stop r^2 f dr.
module;
#include <cassert>
module qchem.Mesh1.Radial.Internal;

LinearRadialMesh::LinearRadialMesh(double start, double stop, int NumPoints)
{
    assert(NumPoints>=2);
    itsR.resize(NumPoints);
    itsW.resize(NumPoints);
    double dr=(stop-start)/(NumPoints-1);
    for (int i=0; i<NumPoints; i++)
    {
        double r=start+i*dr;
        double trap=(i==0 || i==NumPoints-1) ? 0.5*dr : dr;
        itsR[i]=r;
        itsW[i]=trap*r*r;
    }
}
