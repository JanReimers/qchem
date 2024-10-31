// File: AngularIntegrals.H  Handle the angular part of 2-electron ERIs

#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Wigner3j.H"
#include "wignerSymbols/wignerSymbols-cpp.h"
#include "Imp/Misc/IntPower.H"
#include <cassert>

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
//    double w3a=WignerSymbols::wigner3j(la,la,k,0,0,0);
//    double w3c=WignerSymbols::wigner3j(lc,lc,k,0,0,0);
    double w3am=WignerSymbols::wigner3j(la,la,k,ma,-ma,0);
    double w3cm=WignerSymbols::wigner3j(lc,lc,k,mc,-mc,0);
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
//    double w3ab=WignerSymbols::wigner3j(la,lb,k,0,0,0);
    double w3ab_m=WignerSymbols::wigner3j(la,lb,k,ma,-mb,mb-ma);
    return w3ab*w3ab*w3ab_m*w3ab_m;
}

} //namespace 

