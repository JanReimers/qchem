// File: guass.cpp  General gaussian integral.


#include "Imp/Integrals/GaussianIntegrals.H"
#include <iostream>
#include <cassert>

double Oddl (double exp, int n);
double Evenl(double exp, int n);


//##############################################################################
//      /
//  4Pi |  r^l exp(-er^2) r^2 dr
//     /
double GaussianIntegral(double exp, int l)
{
//    return 4*Pi*( l%2 ? Oddl(exp,(l+1)/2) : Evenl(exp,l/2+1) );
    return 4*Pi*( l%2 ? Oddl(exp,(l+1)/2) : Evenl(exp,(l+2)/2 ));
}

double Oddl (double exp, int n)
{
    return n==0 ? 1.0/(2*exp) : Oddl(exp,n-1)*n/exp;
}

double Evenl(double exp, int n)
{
    return n==0 ? sqrt(Pi/exp)/2.0 : Evenl(exp,n-1)*(2*n-1)/(2.0*exp);
//    return n==0 ? sqrt(Pi/exp)/2.0 : Evenl(exp,n-1)*(2*n+1)/(2.0*exp);
}

double GaussianNorm(double e, int l)
{
    double overlap=TwoLPlusOne(l)*GaussianIntegral(2*e,2*l);
//    double overlap=GaussianIntegral(2*e,2*l);
    return 1.0/sqrt(overlap);
}

inline void Swap(double& a, double& b, int& la, int& lb)
{
    double t=a;
    a=b;
    b=t;
    int    it=la;
    la=lb;
    lb=it;
}

//double ExchangeIntegral(double a, double b, int l1, int l2)
//{
//    if (l2 > l1)
//    {
//        int t=l1;
//        l1=l2;
//        l2=t;
//        double td=a;
//        a=b;
//        b=td;
//    }
//    if (l1==0 && l2==0) return GaussianRepulsionIntegral(a,b,l1,l2);
//    double ab=a+b, rab=sqrt(a+b), rab3=ab*rab;
//    if (l1==1 && l2==1) return 3*Pi52/(6*a*b*rab3);
//    double a2=a*a,b2=b*b,rab5=ab*rab3;
//    if (l1==2 && l2==2)
//    {
//        double R0=3*Pi52*(a+ab)*(b+ab)/(2*a2*b2*rab5);
//        double R2=15*Pi52/(2*a*b*rab5);
//        return R0/6+R2/15;
//    }
//    std::cerr << "ExchangeIntegral: Unhandeled exponents in repulsion integral: l1=" << l1 << ", l2=" << l2 << std::endl;
//    return 0;
//}


double GaussianRepulsionIntegral(double ab, double cd, int lab, int lcd)
{
    assert(lab%2==0);
    assert(lcd%2==0);
    if (lab > lcd) Swap(ab,cd,lab,lcd);
    double abcd=ab+cd, r=sqrt(abcd);
    double lambda=2.0*Pi52/(ab*cd*r);

    // <ss|ss>
    if (lab== 0 && lcd== 0) return lambda;
    // <ss|pp>
    if (lab== 0 && lcd== 2) return 3*lambda/6*(1/abcd+2/cd);
    // <pp|pp>
//    double rabcd5=abcd*abcd*r;
//    double ab2=ab*ab, cd2=cd*cd;
//    if (lab== 2 && lcd== 2) return 1.5*Pi52*((ab+2.0*cd)*(cd+2.0*ab))/ (ab2*cd2*rabcd5);
    // We need this for coloumb integrals
    if (lab== 2 && lcd== 2) return 3*lambda/4*(2/(ab*cd)+1/(abcd*abcd));
    // This one is confirmed by summing M&D integrals over {x,y,z}
    // Seems like we need this for exchange integrals.
//    if (lab== 2 && lcd== 2) return 9*lambda/4*(2/(3*ab*cd)+3/(5*abcd*abcd));

// Exchange only
//    double rab3=ab*rab;
//    double a2=a*a, b2=b*b;
//   if (l1==1 && l2==1) return -3.0*Pi52/(6.0*a*b*rab3); //M&D get 3.0*Pi52/(6.0*a2*b2*rab); !?!
// Coulomb only
//    double rab5=ab*rab3;
//    if (l1== 2 && l2== 2) return 1.5*Pi52*((a+2.0*b)*(b+2.0*a))/ (a2*b2*rab5);
//    if (l1== 2 && l2== 2) return lambda/(4.0*a*b)*(2.0/3.0+(3.0*a*b)/(5.0*ab*ab));
//    if (l1== 2 && l2== 0) return Pi52*(3.0*a+2.0*b)/(a2*b *rab3);
//    double a3=a2*a;
//    if (l1== 4 && l2== 0) return 2.00*Pi52/(a*b*rab)*(2.0/a2+1.0/(a*ab)+3.0/4.0/ab2);
//    double ab3=ab2*ab;
//    if (l1== 4 && l2== 2) return 0.75*Pi52*(28*a*b2+8*b2+10*a3+35*a2*b) / (a3*b2*ab3*rab);
//    double ra=sqrt(a), rb=sqrt(b);
//    if (l1== 1 && l2== 1) return 2.00*Pi52*(ra-rab+rb)                  / ab             ;
//    if (l1== 3 && l2== 1) return      Pi52*(rab*(ra+2*rb)-(a+2*b))      / (a2*b*rab)     ;
//  if (l1== 3 && l2== 3) return ;
    std::cerr << "GaussianRepulsionIntegral: Unhandeled exponents in repulsion integral: lab=" << lab << ", lcd=" << lcd << std::endl;
    assert(false);
    return 0;
}

//
//  <ga(r1)gb(r1) 1/r12 gc(r2)gd(r2)>
//

double GaussianRepulsionIntegral(double a, double b, double c, double d, int la, int lb, int lc, int ld)
{
    if (la!=lb || lc!=ld) return 0.0;
    double ab=a+b, cd=c+d;
    if (la==lb && lc==ld) return GaussianRepulsionIntegral(ab,cd,la+lb,lc+ld);
//    if (la>lb) Swap(a,b,la,lb);
//    if (lc>ld) Swap(c,d,lc,ld);
//    double lambda=GaussianRepulsionIntegral(ab,cd,0,0);
//    // Exchange case
//    if ((la==lc && lb==ld) || (la==ld && lb==lc))
//    {
//        //     <ps||ps>      or      <sp||sp> exchange integral
//        if ((la==0 && lb==1) || (la==1 && lb==0))
//        {
//            double abcd=ab+cd;
//            return 3*lambda/(6.0*abcd);
//        }
//    }
//    if (la==lc && lb!=ld) return 0.0;
//    if (la!=lc && lb==ld) return 0.0;

    std::cerr << "GaussianRepulsionIntegral: Unhandeled exponents in repulsion integral: la="
              << la << ", lb=" << lb << ", lc=" << lc << ", ld=" << ld << std::endl;
    assert(false);
    return 0;
}

double GaussianExchangeIntegral(double a, double b, double c, double d, int la, int lb, int lc, int ld)
{
    double ab=a+b, cd=c+d;
    if (la==0 && lb==0 && lc==0 && ld==0) return GaussianRepulsionIntegral(ab,cd,la+lb,lc+ld);
    // TODO maybe we can return zero here.
    if (la==0 && lb==0 && lc==1 && ld==1) return 0.0;//GaussianRepulsionIntegral(ab,cd,la+lb,lc+ld);
    if (la>lb) Swap(a,b,la,lb);
    if (lc>ld) Swap(c,d,lc,ld);
    double lambda=GaussianRepulsionIntegral(ab,cd,0,0);
    double abcd=ab+cd;
    if ((la==lc && lb==ld) || (la==ld && lb==lc))
    {
        //     <ps||ps>      or      <sp||sp> exchange integral
        if ((la==0 && lb==1) || (la==1 && lb==0))
        {
            return 3*lambda/(6.0*abcd);
        }
        //     <pp||pp> exchange integral which is NOT the same as a <pp|pp> Coulomb integral!
        else if ((la==1 && lb==1))
        {
            assert(lc==1 && ld==1);
            return 3*lambda/4*(2/(3*ab*cd)+1/(1*abcd*abcd));
        }
    }
    if (la==lc && lb!=ld) return 0.0;
    if (la!=lc && lb==ld) return 0.0;

    std::cerr << "GaussianExchangeIntegral: Unhandeled exponents in exchange integral: la="
              << la << ", lb=" << lb << ", lc=" << lc << ", ld=" << ld << std::endl;
    assert(false);
    return 0;

}

