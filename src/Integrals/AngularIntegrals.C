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
    return Wigner3j::theW3j(la,k,lb); //What about *(2k+1) ??
}

double Coulomb (int k,int la,int lc,int ma,int mc)
{
    assert(k>=0);
    assert(k<=2*std::min(la,lc));
    int phase=intpow(-1,ma+mc);
    double w3a=WignerSymbols::wigner3j(la,k,la,0,0,0);
    double w3c=WignerSymbols::wigner3j(lc,k,lc,0,0,0);
    double w3am=WignerSymbols::wigner3j(la,k,la,ma,0,-ma);
    double w3cm=WignerSymbols::wigner3j(lc,k,lc,mc,0,-mc);
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
    double w3a=WignerSymbols::wigner3j(la,lb,k,0,0,0);
    double w3b=WignerSymbols::wigner3j(la,lb,k,ma,-mb,mb-ma);
    return w3a*w3a*w3b*w3b;
}

} //namespace 

