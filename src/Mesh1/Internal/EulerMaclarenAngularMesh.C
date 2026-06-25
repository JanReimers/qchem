// File: Internal/EulerMaclarenAngularMesh.C  Euler-Maclaren(theta) x uniform(phi) angular mesh.
// Transplanted VERBATIM from src/Mesh/Internal/EulerMaclarenAngularMesh.C.  m in {1,2,3} sets the
// theta variable transform / clustering.  Weights sum to 4*pi.
module;
#include <cmath>
#include <cassert>
module qchem.Mesh1.Angular.Internal;
import qchem.Math;

namespace qcMesh1
{

EulerMaclarenAngularMesh::EulerMaclarenAngularMesh(int L, int m)
{
    assert(m>=1 && m<=3);
    int numTheta=(L+1)/2;
    int numPhi  =(L+1);
    int numDir  =numTheta*numPhi;
    itsD.resize(numDir);
    itsW.resize(numDir);

    rvec_t thetas(numTheta), Wt(numTheta);
    {
        double del=Pi/numTheta;
        for (int i=0; i<numTheta; i++)
        {
            double q=del*i;
            if (m==1)
            {
                thetas[i]=q;
                Wt[i]    =del*std::sin(thetas[i]);
            }
            if (m==2)
            {
                thetas[i]=q*q*(3*Pi-2*q)/(Pi*Pi);
                Wt[i]    =del*std::sin(thetas[i])*6.0/(Pi*Pi)*q*(Pi-q);
            }
            if (m==3)
            {
                thetas[i]=q*q*q*(10*Pi*Pi-15*Pi*q+6*q*q)/(Pi*Pi*Pi*Pi);
                Wt[i]    =del*std::sin(thetas[i])*30.0/(Pi*Pi*Pi*Pi)*q*q*(Pi-q)*(Pi-q);
            }
        }
    }
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
