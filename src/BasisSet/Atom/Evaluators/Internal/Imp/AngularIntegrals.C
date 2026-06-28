// File: BasisSet/Atom/Evaluators/Internal/Imp/AngularIntegrals.C  Handle the angular part of 2-electron ERIs
module;
#include <cassert>
#include <iostream>
module qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;
import qchem.BasisSet.Atom.Evaluators.Internal.Wigner3j;
import qchem.Math;

using std::cout;
using std::endl;

namespace AngularIntegrals
{

double Exchange(int k,int la,int lb)
{
    assert(k>=std::abs(la-lb));
    assert(k<=la+lb);
    return FourPi2*Wigner3j::w3j.sq(la,lb,k);   // (3j)^2 directly (no sqrt)
}

double Direct  (int k,int la,int lc,int ma,int mc)
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
    // (3j)^2 * (3j_m)^2 directly (no sqrt) -- exchange squares the symbols.
    return FourPi2*Wigner3j::w3j.sq(la,lb,k)*Wigner3j::w3j.sq(la,lb,k,ma,-mb);
}

rvec11_t Direct (int la,int lb)
{    
    rvec11_t Ak(0.0);
    Ak[0]=FourPi2;
    return Ak;
}


rvec11_t Exchange(int la,int lb)
{    
    assert(la>=0);
    assert(lb>=0);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    rvec11_t Ak(0.0);
    assert(kmax<=Ak.size());
    for (int k=kmin;k<=kmax;k+=2)
    {
        assert((k+la+lb)%2==0);
        Ak[k]=Exchange(k,la,lb); //What about *(2k+1) ??
    }
    return Ak;
}


rvec11_t Direct  (int la,int lc,int ma,int mc)
{    
    rvec11_t Ak(0.0);
    int kmax=2*std::min(la,lc);
    assert(kmax<=Ak.size());
    int phase=intpow(-1,ma+mc);
    for (int k=0;k<=kmax;k+=2)
    {
        double w3a=Wigner3j::w3j(la,la,k);
        double w3c=Wigner3j::w3j(lc,lc,k);
        double w3am=Wigner3j::w3j(la,la,k,ma,-ma);
        double w3cm=Wigner3j::w3j(lc,lc,k,mc,-mc);
        Ak[k]= FourPi2*phase*(2*la+1)*(2*lc+1)*w3a*w3am*w3c*w3cm;
    }
    return Ak;
}

rvec11_t Exchange(int la,int lb,int ma,int mb)
{    
    assert(la>=0);
    assert(lb>=0);
    assert(ma>=-la);
    assert(ma<= la);
    assert(mb>=-lb);
    assert(mb<= lb);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    rvec11_t Ak(0.0);
    assert(kmax<=Ak.size());
    for (int k=kmin;k<=kmax;k+=2)
    {
        assert((k+la+lb)%2==0);
        // (3j)^2 * (3j_m)^2 directly (no sqrt) -- exchange squares the symbols.
        Ak[k]=FourPi2*(2*la+1)*(2*lb+1)*Wigner3j::w3j.sq(la,lb,k)*Wigner3j::w3j.sq(la,lb,k,ma,-mb); //What about *(2k+1) ??
    }
    return Ak;
}

} //namespace 

