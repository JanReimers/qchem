// File: Internal/GaussLegendreAngularMesh.C  Gauss-Legendre(theta) x uniform(phi) angular mesh.
// Transplanted from src/Mesh/Internal/GaussLegendreAngularMesh.C, now using the SHARED
// qchem.Mesh1.GaussLegendre (whose indexing bug is fixed) instead of a private GL copy.
// Weights: sum_theta Wt ~ 2, sum_phi Wp = 2*pi, so sum = 4*pi.
module;
#include <cmath>
#include <cassert>
module qchem.Mesh1.Angular.Internal;
import qchem.Mesh1.GaussLegendre;
import qchem.Math;

namespace qcMesh1
{

GaussLegendreAngularMesh::GaussLegendreAngularMesh(int L)
{
    int numTheta=(L+1)/2;
    int numPhi  =(L+1);
    int numDir  =numTheta*numPhi;
    itsD.resize(numDir);
    itsW.resize(numDir);

    // theta nodes from Gauss-Legendre on cos(theta) in [-1,1].
    GaussLegendre gl(numTheta,-1.0,1.0);
    rvec_t thetas(numTheta), Wt(numTheta);
    for (int i=0; i<numTheta; i++)
    {
        thetas[i]=std::acos(gl.x[i]);
        Wt[i]    =gl.w[i];
    }
    // phi nodes uniform on [0,2*pi).
    double delPhi=2*Pi/numPhi;

    int k=0;
    for (int it=0; it<numTheta; it++)
        for (int ip=0; ip<numPhi; ip++,k++)
        {
            double t=thetas[it], p=delPhi*ip;
            itsD[k]=rvec3_t(std::sin(t)*std::sin(p), std::sin(t)*std::cos(p), std::cos(t));
            itsW[k]=Wt[it]*delPhi;
        }
    assert(k==numDir);
}

} //namespace qcMesh1
