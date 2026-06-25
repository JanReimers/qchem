// File: Internal/MHLRadialMesh.C  Murray-Handy-Laming radial mesh builder.  Transplanted VERBATIM
// from src/Mesh/Internal/MHLRadialMesh.C.
//   r_i = alpha (x/(1-x))^m,   w_i = (del alpha m x^{m-1}/(1-x)^{m+1}) r_i^2,   x = i/N.
// The r^2 jacobian is folded into w_i.
module;
#include <utility>
module qchem.Mesh1.Radial;
import qchem.Math;

namespace qcMesh1
{

RadialMesh MHLRadial(int NumPoints, int m, double alpha)
{
    rvec_t R(NumPoints), W(NumPoints);
    double del=1.0/NumPoints;
    for (int i=0; i<NumPoints; i++)
    {
        double x=i*del;
        double r=alpha*intpow(x,m)/intpow(1.0-x,m);
        R[i]=r;
        W[i]=del*alpha*m*intpow(x,m-1)/intpow(1.0-x,m+1)*r*r;
    }
    return RadialMesh(std::move(R), std::move(W));
}

} //namespace qcMesh1
