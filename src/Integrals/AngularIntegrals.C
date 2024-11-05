// File: AngularIntegrals.H  Handle the angular part of 2-electron ERIs

#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector.h"
#include <cassert>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

namespace AngularIntegrals
{
    
double Exchange(int k,int la,int lb)
{
    assert(k>=std::abs(la-lb));
    assert(k<=la+lb);
    double wabk=Wigner3j::w3j(la,lb,k);
    return wabk*wabk; 
}

double Coulomb (int k,int la,int lc,int ma,int mc)
{
    assert(k>=0);
    assert(k<=2*std::min(la,lc));
    int phase=intpow(-1,ma+mc);
    double w3a=Wigner3j::w3j(la,la,k);
    double w3c=Wigner3j::w3j(lc,lc,k);
    double w3am=Wigner3j::w3j(la,la,k,ma,-ma);
    double w3cm=Wigner3j::w3j(lc,lc,k,mc,-mc);
//    if ((la==2 || lc==2) && (fabs(w3a) != fabs(w3am)))
//        cout << "la,lc,k,w3a,w3c,w3am,w3cm=" << la << " " << lc << " " << k << " " <<  w3a << " " << w3c << " " << w3am << " " << w3cm << endl;

    return phase*w3a*w3am*w3c*w3cm;
}

double Exchange(int k,int la,int lb,int ma,int mb)
{
    assert(la>=0);
    assert(lb>=0);
    assert(ma>=-la);
    assert(ma<= la);
    assert(mb>=-lb);
    assert(mb<= lb);
    assert(k>=std::abs(la-lb));
    assert(k<=la+lb);
    double w3ab=Wigner3j::w3j(la,lb,k);
    double w3ab_m=Wigner3j::w3j(la,lb,k,ma,-mb);
//    if ((la==2 || lb==2))
//        cout << "la,lb,k,w3ab,w3abm = (" << la << " " << lb << " " << k << ")    " <<  w3ab << "    " <<  w3ab_m << endl;
    return w3ab*w3ab*w3ab_m*w3ab_m;
}

RVec Coulomb (int la,int lc,int ma,int mc)
{    
    RVec Ak(la+lc+1);
    int kmax=2*std::min(la,lc);
    int phase=intpow(-1,ma+mc);
    int i=1;
    for (int k=0;k<=kmax;k+=2,i++)
    {
        double w3a=Wigner3j::w3j(la,la,k);
        double w3c=Wigner3j::w3j(lc,lc,k);
        double w3am=Wigner3j::w3j(la,la,k,ma,-ma);
        double w3cm=Wigner3j::w3j(lc,lc,k,mc,-mc);
        Ak(i)= phase*w3a*w3am*w3c*w3cm;
    }
    return Ak;
}

} //namespace 

