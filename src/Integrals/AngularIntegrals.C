// File: AngularIntegrals.H  Handle the angular part of 2-electron ERIs

#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Misc/IntPower.H"
#include <cassert>
#include <cmath>

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
    return w3ab*w3ab*w3ab_m*w3ab_m;
}

} //namespace 

