// File:  AuxillaryFJ.C  Class for calculating the auxiliary function.

#include "PolarizedGaussian/AuxillaryFJ.H"
#include "Common/IntPower.H"
#include "oml/vector.h"
#include <cmath>
#include <cassert>

import Common.Constants;


double Factorial[]= {1,1,2,6,24,120,720};

void AuxillaryFJ::GetFjAt(double T, Vector<double>& Fj) const
{
    assert(Fj.GetLimits().Low ==0        ); //Illegal lower limit for Fj.
    assert(Fj.GetLimits().High<=thejMax+6); //Lookup table not large enough.
    assert(T>=0                          );

    Vector<double>::Subscriptor s(Fj);
    int jmax=Fj.GetLimits().High;
    double et=exp(-T);

    if(T<12)
    {
        int it=(int)(T*10.0+0.5);
        double Tstar=(double)it*0.1;
        double dt=Tstar-T;
        double* FT=&theLookUp[it][0];
        double fjmax=FT[jmax];
        for (int k=1; k<=6; k++) fjmax+=FT[k+jmax]*uintpow(dt,k)/Factorial[k];
        s(jmax)=fjmax;
        for (int j=jmax-1; j>=0; j--) s(j)=(2*T*s(j+1)+et)/(2*j+1);
    }
    if(T>=12 && T < 30.0)
    {
        double g=0;
        if(         T<15) g = 0.4999489092 - 0.2473631686/T + 0.321180909/(T*T) - 0.3811559346/(T*T*T);
        if(T>=15 && T<18) g = 0.4998436875 - 0.24249438  /T + 0.24642845 /(T*T);
        if(T>=18 && T<24) g = 0.499093162  - 0.2152832   /T;
        if(T>=24        ) g = 0.49;
        assert(g!=0.0); // g never got assigned.
        s(0) = 0.5*sqrt(Pi/T) - et*g/T;
        for (int j=0; j<jmax; j++) s(j+1)=((2*j+1)*s(j)-et)/(2*T);
    }
    if (T>=30.0)
    {
        s(0)=0.5*sqrt(Pi/T);
        for (int j=0; j<jmax; j++) s(j+1)=((2*j+1)*s(j)-et)/(2*T);
    }
}
