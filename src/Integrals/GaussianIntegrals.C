// File: guass.cpp  General gaussian integral.


#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/GaussianRadialIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Factorials.H"

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
    return 1.0/sqrt(GaussianIntegral(2*e,2*l));
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

double GaussianRadialIntegrals::FourPi2=4*4*Pi*Pi;


GaussianRadialIntegrals::GaussianRadialIntegrals(double e_ab, double e_cd)
    : a(e_ab)
    , c(e_cd)
{
}

double GaussianRadialIntegrals::R(int k,int la, int lb, int lc, int ld) const
{
    int Lab_p=la+lb+k;
    int Lab_m=la+lb-k;
    int Lcd_p=lc+ld+k;
    int Lcd_m=lc+ld-k;
    //Check that everything is even.
    assert(Lab_p%2==0);
    assert(Lab_m%2==0);
    assert(Lcd_p%2==0);
    assert(Lcd_m%2==0);
    assert(Lab_m>=0);
    assert(Lcd_m>=0);

    assert(Lab_p+1<=qchem::NMax);
    assert(Lcd_p+1<=qchem::NMax);
    double afact=qchem::DFact[Lcd_p+1]/qchem::Twon[Lcd_p/2]; //These ab and cd are reversed on purpose.
    double cfact=qchem::DFact[Lab_p+1]/qchem::Twon[Lab_p/2];
    double Iab=Dab(Lab_m/2,Lcd_p+3);
    double Icd=Dcd(Lcd_m/2,Lab_p+3);
    //if(Lab_p>0 || Lcd_p>0 )
    {
        //std::cout <<  "L+=" << Lab_p << ", " << Lcd_p << ", L-=" << Lab_m << ", " << Lcd_m << std::endl;
        //std::cout << k << "  " <<  afact << " * " << Iab << "   +   " << cfact << " * " << Icd << std::endl;
    }
    return Pi12/8*(afact*Iab+cfact*Icd);
}
double GaussianRadialIntegrals::Coulomb(int la, int lb, int lc, int ld) const
{
    assert(la==lb);
    assert(lc==ld);
    return FourPi2*R(0,la,lb,lc,ld);
}

double GaussianRadialIntegrals::DoExchangeSum(int la, int lb, int lc, int ld) const
{
    assert(la==lc);
    assert(lb==ld);
    assert(la>=0);
    assert(lb>=0);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    double ret=0.0;
    for (int k=kmin;k<=kmax;k+=2)
    {
        ret+=R(k,la,lb,la,lb)*AngularIntegrals::Exchange(k,la,lb);
    }
    return FourPi2*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
}


double GaussianRadialIntegrals::D(double _a, int m,int n) const
{
    double ret=0.0;
    switch (m)
    {
    case 0:
    {
        ret=I(_a,0,n);
        break;
    }
    case 1:
    {
        ret=I(_a,0,n)*I(_a,1,n);
        break;
    }
    case 2:
    {
        double I1=I(_a,1,n);
        ret=I(_a,0,n)*(I1*I1+I(_a,2,n));
        break;
    }
    case 3:
    {
        double I1=I(_a,1,n);
        double I2=I(_a,2,n);
        ret=I(_a,0,n)*(I1*I1*I1+3*I1*I2+I(_a,3,n));
        break;
    }
    default :
    {
        assert(false);
    }
    }
    return ret;
}

double GaussianRadialIntegrals::I(double _a, int m,int n) const
{
    assert(n>0);
    assert(m>=0);
    assert(m<=qchem::NMax);
    assert(_a==a || _a==c);
    if (m==0)
    {
        return 1.0/(_a*sqrt(pow(a+c,n)));
    }
//    else if (m==1)
//        return n/(2*(a+c))+1/_a;
//    else if (m==2)
//        return n/(2*(a+c)*(a+c))+1/(_a*_a);
    else
        return qchem::Fact[m-1]*(n/(2*pow((a+c),m))+1/pow(_a,m));
}

double GaussianRadialIntegrals::Dab(int m,int n) const
{
    return D(a,m,n);
}

double GaussianRadialIntegrals::Dcd(int m,int n) const
{
    return D(c,m,n);
}


