// File: Internal/GaussLegendreAngularMesh.C  Gauss-Legendre(theta) x uniform(phi) angular builder.
// Uses the SHARED qchem.Mesh1.GaussLegendre (whose indexing bug is fixed).
// Weights: sum_theta Wt ~ 2, sum_phi Wp = 2*pi, so sum = 4*pi.
module;
#include <cmath>
#include <cassert>
#include <utility>
module qchem.Mesh1.Angular;
import qchem.Mesh1.GaussLegendre;
import qchem.Math;

namespace qcMesh1
{

AngularMesh GaussLegendreAngular(int L)
{
    int numTheta=(L+1)/2;
    int numPhi  =(L+1);
    int numDir  =numTheta*numPhi;
    rvec3vec_t D(numDir);
    rvec_t     W(numDir);

    // theta nodes from Gauss-Legendre on cos(theta) in [-1,1].
    GaussLegendre gl(numTheta,-1.0,1.0);
    rvec_t thetas(numTheta), Wt(numTheta);
    for (int i=0; i<numTheta; i++)
    {
        thetas[i]=std::acos(gl.x[i]);
        Wt[i]    =gl.w[i];
    }
    double delPhi=2*Pi/numPhi;            // phi nodes uniform on [0,2*pi)

    int k=0;
    for (int it=0; it<numTheta; it++)
        for (int ip=0; ip<numPhi; ip++,k++)
        {
            double t=thetas[it], p=delPhi*ip;
            D[k]=rvec3_t(std::sin(t)*std::sin(p), std::sin(t)*std::cos(p), std::cos(t));
            W[k]=Wt[it]*delPhi;
        }
    assert(k==numDir);
    return AngularMesh(std::move(D), std::move(W));
}

} //namespace qcMesh1
