// File: Internal/LogRadialMesh.C  Logarithmic radial mesh.  Transplanted VERBATIM from
// src/Mesh/Internal/LogRadialMesh.C (same node/weight formulae, same special first/last weights:
// the first node integrates the whole inner sphere, the last only a half annulus).
module;
module qchem.Mesh1.Radial.Internal;
import qchem.Math;

namespace qcMesh1
{

LogRadialMesh::LogRadialMesh(double start, double stop, int NumPoints)
{
    itsR.resize(NumPoints);
    itsW.resize(NumPoints);
    double q  = exp((log(stop)-log(start))/(NumPoints-1));
    double sq = sqrt(q);
    double wq = 1.0/3.0 * ( Cube(sq) - Cube(1.0/sq) );
    double r  = start;
    itsR[0]=r;
    itsW[0]=1.0/3.0 * Cube(start*sq);            // whole sphere instead of an annulus
    r*=q;
    for (int i=1; i<NumPoints-1; i++)
    {
        itsR[i]=r;
        itsW[i]=Cube(r)*wq;
        r*=q;
    }
    itsR[NumPoints-1]=r;                          // == stop
    itsW[NumPoints-1]=1.0/3.0 * (Cube(stop)-Cube(stop/sq)); // only half annulus
}

} //namespace qcMesh1
