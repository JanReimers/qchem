//  File: LogRadialMesh.C  Logarithmic radial mesh implementation



#include "LogRadialMesh.H"
#include <iostream>
#include <cmath>

import Common.Constants;

LogRadialMesh::LogRadialMesh(double start, double stop, index_t NumPoints)
{
    double q   = exp((log(stop)-log(start))/(NumPoints-1));
    double sq  = sqrt(q);
    double wq  = 1.0/3.0 * ( Cube(sq) - Cube(1.0/sq) );
    double r   = start;
    push_back(r,1.0/3.0 *  Cube(start*sq)); //Do whole sphere instead of anulus.
    r*=q;
    for(int i=2;i<NumPoints;i++)
    {
        push_back(r,Cube(r)*wq);
        r*=q;
    }
    push_back(r,1.0/3.0 * (Cube(stop)-Cube(stop/sq))); //Do only half anulus.
//    Initialize(R,W);
}

