//  File: GaussLegendreAngularMesh.C  GaussLegendre style angular mesh implementation



#include "Mesh/AngularMesh/GaussLegendreAngularMesh.H"
#include "Misc/DFTDefines.H"
#include "oml/vector.h"
#include <iostream>
#include <cmath>
#include <cassert>

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
