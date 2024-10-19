//  File: MHLRadialMesh.C  MHLarithmic radial mesh implementation



#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Misc/DFTDefines.H"
#include "Misc/IntPower.H"
#include <iostream>
#include <cmath>

MHLRadialMesh::MHLRadialMesh(index_t NumPoints, int m, double alpha)
    : RadialMeshImplementation()
{
    Vector<double> R(NumPoints);
    Vector<double> W(NumPoints);

    Vector<double>::iterator ir(R.begin());
    Vector<double>::iterator iw(W.begin());
    double del=1.0/NumPoints;
    for(index_t i=0; i<NumPoints; i++,ir++,iw++)
    {
        double x=i*del;
        *ir=alpha*intpow(x,m)/intpow(1.0-x,m);
        *iw=del*alpha*m*intpow(x,m-1)/intpow(1.0-x,m+1)*(*ir)*(*ir);
    }

#if DEBUG_OUTPUT
    double rmax=R(NumPoints);
    cout << "Sum of weigths/Vol(rmax) = " << Sum(W)/(4.0/3.0*Pi*rmax*rmax*rmax) << std::endl;
#endif

    Initialize(R,W);
}

