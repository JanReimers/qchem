// File: AngularIntegrals.C  Handle the angular part of 2-electron ERIs

#include "AngularIntegrals.H"
#include "Wigner3j.H"
#include "Common/IntPower.H"
#include "oml/vector.h"
#include <cassert>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

namespace AngularIntegrals
{

double FourPi2=4*4*pi*pi;
    
double Exchange(int k,int la,int lb)
{
    assert(k>=std::abs(la-lb));
    assert(k<=la+lb);
    double wabk=Wigner3j::w3j(la,lb,k);
    return FourPi2*wabk*wabk; 
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

    return FourPi2*phase*w3a*w3am*w3c*w3cm;
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
    return FourPi2*w3ab*w3ab*w3ab_m*w3ab_m;
}

RVec Coulomb(int la,int lb)
{    
    RVec Ak(1);
    Ak(1)=FourPi2;
    return Ak;
}


RVec Exchange(int la,int lb)
{    
    assert(la>=0);
    assert(lb>=0);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    int N=(kmax-kmin)/2+1;
    RVec Ak(N);
    int i=1;
    for (int k=kmin;k<=kmax;k+=2)
    {
        assert((k+la+lb)%2==0);
        Ak(i++)=Exchange(k,la,lb); //What about *(2k+1) ??
    }
    return Ak;
}


RVec Coulomb (int la,int lc,int ma,int mc)
{    
    RVec Ak(la+lc+1,0.0);
    int kmax=2*std::min(la,lc);
    int phase=intpow(-1,ma+mc);
    int i=1;
    for (int k=0;k<=kmax;k+=2,i++)
    {
        double w3a=Wigner3j::w3j(la,la,k);
        double w3c=Wigner3j::w3j(lc,lc,k);
        double w3am=Wigner3j::w3j(la,la,k,ma,-ma);
        double w3cm=Wigner3j::w3j(lc,lc,k,mc,-mc);
        Ak(i)= FourPi2*phase*(2*la+1)*(2*lc+1)*w3a*w3am*w3c*w3cm;
    }
    return Ak;
}

RVec Exchange(int la,int lb,int ma,int mb)
{    
    assert(la>=0);
    assert(lb>=0);
    assert(ma>=-la);
    assert(ma<= la);
    assert(mb>=-lb);
    assert(mb<= lb);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    int N=(kmax-kmin)/2+1;
    RVec Ak(N,0.0);
    int i=1;
    for (int k=kmin;k<=kmax;k+=2)
    {
        assert((k+la+lb)%2==0);
        double w3ab=Wigner3j::w3j(la,lb,k);
        double w3ab_m=Wigner3j::w3j(la,lb,k,ma,-mb);
        Ak(i++)=FourPi2*(2*la+1)*(2*lb+1)*w3ab*w3ab*w3ab_m*w3ab_m; //What about *(2k+1) ??
    }
    return Ak;
}

} //namespace 

