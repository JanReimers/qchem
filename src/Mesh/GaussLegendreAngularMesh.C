//  File: GaussLegendreAngularMesh.C  GaussLegendre style angular mesh implementation
module;
#include "oml/vector.h"
#include <iostream>
#include <cmath>
#include <cassert>

export module Mesh.GaussLegendreAngularMesh;
export import Mesh;

import Common.Constants;

export class GaussLegendreAngularMesh :  public  Mesh
{
public:
    GaussLegendreAngularMesh(int L, int m);
};


void GaussLegendre(double x1, double x2, Vector<double>& x, Vector<double>& w, int n);

GaussLegendreAngularMesh::GaussLegendreAngularMesh(int L, int)
{
    int numTheta = (L+1)/2;
    int numPhi   = (L+1);
    int numDir   = numTheta*numPhi;

    Vector<RVec3>  D(numDir);
    Vector<double> W(numDir);
//
//  Get a bunch of thetas with weights.
//
    Vector<double> thetas(numTheta);
    Vector<double> Wt(numTheta);
    GaussLegendre(-1.0,1.0,thetas,Wt,numTheta);
    thetas=acos(thetas);
#if DEBUG_OUTPUT
    cout << "Sum of theta weigths = " << Sum(Wt) << std::endl;
#endif
//
//  Get a bunch of phi's with weights.
//
    Vector<double> phis(numPhi);
    Vector<double> Wp  (numPhi);
    {
        Vector<double>::iterator it(phis.begin());
        Vector<double>::iterator iw(Wp.begin());
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
    Vector<RVec3 >::iterator d(D.begin());
    Vector<double>::iterator w(W.begin());

    Vector<double>::const_iterator t (thetas.begin());
    Vector<double>::const_iterator wt(Wt.begin());
    for(; t; t++,wt++)
    {
        Vector<double>::const_iterator p (phis.begin());
        Vector<double>::const_iterator wp(Wp.begin());
        for (; p; p++,wp++,d++,w++)
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
    for (auto i:D.indices()) push_back(D(i),W(i));
}

void GaussLegendre(double x1, double x2, Vector<double>& x, Vector<double>& w, int n)
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
        x(i    )=xm-xl*z;
        x(n+1-i)=xm+xl*z;
        w(i)=2.0*xl/((1.0-z*z)*pp*pp);
        w(n+1-i)=w(i);
    }
}
