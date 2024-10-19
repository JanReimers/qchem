//  File: EulerMaclarenAngularMesh.C  EulerMaclaren style angular mesh implementation



#include "Imp/Mesh/EulerMaclarenAngularMesh.H"
#include "Imp/Misc/DFTDefines.H"
#include "oml/vector.h"
#include <iostream>
#include <cmath>
#include <cassert>

EulerMaclarenAngularMesh::EulerMaclarenAngularMesh(int L, int m) 
{
    assert(m>=1 && m<=3);
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
    {
        Vector<double>::iterator it(thetas.begin());
        Vector<double>::iterator iw(Wt.begin());
        double del=Pi/numTheta;
        for (int i=0; i<numTheta; i++,it++,iw++)
        {
            double q=del*i;
            if (m==1)
            {
                *it=q;
                *iw=del*sin(*it);
            }
            if (m==2)
            {
                *it=q*q*(3*Pi-2*q)/(Pi*Pi);
                *iw=del*sin(*it)*6.0/(Pi*Pi)*q*(Pi-q);
            }
            if (m==3)
            {
                *it=q*q*q*(10*Pi*Pi-15*Pi*q+6*q*q)/(Pi*Pi*Pi*Pi);
                *iw=del*sin(*it)*30.0/(Pi*Pi*Pi*Pi)*q*q*(Pi-q)*(Pi-q);
            }
        }
    }
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
#if DEBUG_OUTPUT
    cout << "Sum of phi weigths = " << Sum(Wp) << std::endl;
#endif
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

