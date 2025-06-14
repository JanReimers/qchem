//  File: MHLRadialMesh.C  MHLarithmic radial mesh implementation



#include "Imp/Mesh/MHLRadialMesh.H"
#include "Common/DFTDefines.H"
#include "Common/IntPower.H"
#include <iostream>
#include <cmath>

MHLRadialMesh::MHLRadialMesh(index_t NumPoints, int m, double alpha)
{
    double del=1.0/NumPoints;
    for(index_t i=0; i<NumPoints; i++)
    {
        double x=i*del;
        double r=alpha*intpow(x,m)/intpow(1.0-x,m);
        double w=del*alpha*m*intpow(x,m-1)/intpow(1.0-x,m+1)*r*r;
        push_back(r,w);
    }
}

