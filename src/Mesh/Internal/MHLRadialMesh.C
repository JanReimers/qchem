//  File: MHLRadialMesh.C  MHLarithmic radial mesh implementation
module;
#include <cmath>

module qchem.Mesh.Internal.RadialTypes;
import Common.IntPow;
import qchem.RadialMesh;
// import oml;



MHLRadialMesh::MHLRadialMesh(int NumPoints, int m, double alpha)
{
    double del=1.0/NumPoints;
    for(int i=0; i<NumPoints; i++)
    {
        double x=i*del;
        double r=alpha*intpow(x,m)/intpow(1.0-x,m);
        double w=del*alpha*m*intpow(x,m-1)/intpow(1.0-x,m+1)*r*r;
        push_back(r,w);
    }
}

