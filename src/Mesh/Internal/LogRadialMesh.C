// File: Internal/LogRadialMesh.C  Logarithmic radial mesh builder.  Transplanted VERBATIM from
// src/Mesh/Internal/LogRadialMesh.C (same node/weight formulae, same special first/last weights:
// the first node integrates the whole inner sphere, the last only a half annulus).
module;
#include <utility>
module qchem.Mesh.Radial;
import qchem.Math;

namespace qcMesh
{

RadialMesh LogRadial(double start, double stop, int NumPoints)
{
    rvec_t R(NumPoints), W(NumPoints);
    double q  = exp((log(stop)-log(start))/(NumPoints-1));
    double sq = sqrt(q);
    double wq = 1.0/3.0 * ( Cube(sq) - Cube(1.0/sq) );
    double r  = start;
    R[0]=r;
    W[0]=1.0/3.0 * Cube(start*sq);            // whole sphere instead of an annulus
    r*=q;
    for (int i=1; i<NumPoints-1; i++)
    {
        R[i]=r;
        W[i]=Cube(r)*wq;
        r*=q;
    }
    R[NumPoints-1]=r;                          // == stop
    W[NumPoints-1]=1.0/3.0 * (Cube(stop)-Cube(stop/sq)); // only half annulus
    return RadialMesh(std::move(R), std::move(W));
}

} //namespace qcMesh
