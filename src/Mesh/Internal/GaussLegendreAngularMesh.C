//  File: GaussLegendreAngularMesh.C  GaussLegendre style angular mesh implementation
module;
#include <iostream>
#include <cmath>
#include <cassert>
#include <blaze/Math.h>

module qchem.Mesh.Internal.Types;
import Common.Constants;

void GaussLegendre(double x1, double x2, rvec_t& x, rvec_t& w, int n)
{
    const double eps=3.0E-14;
    double pp=0, z=0, z1=0;
    int m=(n+1)/2;
    double xm=0.5*(x2+x1), xl=0.5*(x2-x1);
    for (int i=1; i<=m; i++)
    {
        z=cos(Pi*(i-0.25)/(n+0.5));
        do
        {
            double p1=1.0,p2=0.0,p3=0.0;
            for (int j=1; j<=n; j++)
            {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        }
        while (fabs(z-z1) > eps);
        x[i-1  ]=xm-xl*z;
        x[n+1-i+1]=xm+xl*z;
        w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i+1]=w[i-1];
    }
}

GaussLegendreAngularMesh::GaussLegendreAngularMesh(int L, int)
{
    int numTheta = (L+1)/2;
    int numPhi   = (L+1);
    int numDir   = numTheta*numPhi;

    rvec3vec_t D(numDir);
    rvec_t     W(numDir);
//
//  Get a bunch of thetas with weights.
//
    rvec_t thetas(numTheta);
    rvec_t Wt(numTheta);
    GaussLegendre(-1.0,1.0,thetas,Wt,numTheta);
    thetas=blaze::acos(thetas);
#if DEBUG_OUTPUT
    cout << "Sum of theta weigths = " << Sum(Wt) << std::endl;
#endif
//
//  Get a bunch of phi's with weights.
//
    rvec_t phis(numPhi);
    rvec_t Wp  (numPhi);
    {
        auto it(phis.begin());
        auto iw(Wp.begin());
        double del=2*Pi /numPhi;
        for (int i=0; i<numPhi; i++,it++,iw++)
        {
            double q=del*i;
            *it=q;
            *iw=del;
        }
    }
//  cout << "Sum of phi weigths = " << Sum(Wp) << std::endl;
//
//  Now take the direct product of the
//
    auto d(D.begin());
    auto w(W.begin());

    auto t (thetas.begin());
    auto wt(Wt.begin());
    for(; t!=thetas.end(); t++,wt++)
    {
        auto p (phis.begin());
        auto wp(Wp.begin());
        for (; p!=phis.end(); p++,wp++,d++,w++)
        {
            (*d).x=sin(*t)*sin(*p);
            (*d).y=sin(*t)*cos(*p);
            (*d).z=cos(*t);
            *w    =(*wt)*(*wp);
        }
    }
#if DEBUG_OUTPUT
    cout << "Sum of weigths/4Pi = " << Sum(W)/4/Pi << std::endl;
#endif
    assert(W.size()==D     .size());
    for (auto i:iv_t(0,D.size())) push_back(D[i],W[i]);
}

