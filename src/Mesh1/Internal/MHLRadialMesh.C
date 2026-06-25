// File: Internal/MHLRadialMesh.C  Murray-Handy-Laming radial mesh.  Transplanted VERBATIM
// from src/Mesh/Internal/MHLRadialMesh.C.
//   r_i = alpha (x/(1-x))^m,   w_i = (del alpha m x^{m-1}/(1-x)^{m+1}) r_i^2,   x = i/N.
// The r^2 jacobian is folded into w_i.
module;
module qchem.Mesh1.Radial.Internal;
import qchem.Math;

MHLRadialMesh::MHLRadialMesh(int NumPoints, int m, double alpha)
{
    itsR.resize(NumPoints);
    itsW.resize(NumPoints);
    double del=1.0/NumPoints;
    for (int i=0; i<NumPoints; i++)
    {
        double x=i*del;
        double r=alpha*intpow(x,m)/intpow(1.0-x,m);
        double w=del*alpha*m*intpow(x,m-1)/intpow(1.0-x,m+1)*r*r;
        itsR[i]=r;
        itsW[i]=w;
    }
}
