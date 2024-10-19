//  File: LogRadialMesh.C  Logarithmic radial mesh implementation



#include "Mesh/RadialMesh/LogRadialMesh.H"
#include "Misc/DFTDefines.H"
#include <iostream>
#include <cmath>

LogRadialMesh::LogRadialMesh(double start, double stop, index_t NumPoints)
    : RadialMeshImplementation()
{
    Vector<double> R(NumPoints);
    Vector<double> W(NumPoints);
    double q   = exp((log(stop)-log(start))/(NumPoints-1));
    double sq  = sqrt(q);
    double wq  = 1.0/3.0 * ( Cube(sq) - Cube(1.0/sq) );
    double r   = start;

    Vector<double>::iterator ir(R.begin());
    Vector<double>::iterator iw(W.begin());
    for(;ir!=R.end()&&iw!=W.end(); ir++,iw++)
    {
        *ir=r;
        *iw=Cube(r)*wq;
        r*=q;
    }
    W(1)         = 1.0/3.0 *  Cube(start*sq); //Do whole sphere instead of anulus.
    W(NumPoints) = 1.0/3.0 * (Cube(stop)-Cube(stop/sq)); //Do only half anulus.

#if DEBUG_OUTPUT
    double rmax=R(NumPoints);
    cout << "Sum of weigths/Vol(rmax) = " << Sum(W)/(4.0/3.0*Pi*rmax*rmax*rmax) << std::endl;
#endif

    Initialize(R,W);
}

