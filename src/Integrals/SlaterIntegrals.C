#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Wigner3j.H"
#include "Misc/DFTDefines.H"
#include <cassert>
#include <cmath>
#include <iostream>

const int NMax=13; //Should enable up to f orbitals.

double SlaterIntegrals::DFact[NMax+1]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
double SlaterIntegrals::Fact[NMax+1]; //factorials lookup table.
double SlaterIntegrals::Twon[NMax+1];  //2^n lookup table.

SlaterIntegrals::SlaterIntegrals(double e_ab, double e_cd)
    : a(e_ab)
    , c(e_cd)
{
    if (DFact[0]!=1.0)
    {
        DFact[0]=1.0;
        DFact[1]=1.0;
        Fact[0]=1.0;
        Twon[0]=1.0;
        Fact[1]=1.0;
        Twon[1]=2.0;
        for (int n=2; n<=NMax; n++)
        {
            DFact[n]=DFact[n-2]*n;
            Fact[n]=Fact[n-1]*n;
            Twon[n]=Twon[n-1]*2;
        }
    }
}

double SlaterIntegrals::operator()(int l,int la, int lb, int lc, int ld) const
{
    int Lab_p=la+lb+l;
    int Lab_m=la+lb-l;
    int Lcd_p=lc+ld+l;
    int Lcd_m=lc+ld-l;
    //Check that everything is even.
    assert(Lab_p%2==0);
    assert(Lab_m%2==0);
    assert(Lcd_p%2==0);
    assert(Lcd_m%2==0);
    assert(Lab_m>=0);
    assert(Lcd_m>=0);

    assert(Lab_p+1<=NMax);
    assert(Lcd_p+1<=NMax);
    double afact=DFact[Lcd_p+1]/Twon[Lcd_p/2]; //These ab and cd are reversed on purpose.
    double cfact=DFact[Lab_p+1]/Twon[Lab_p/2];
    double Iab=Dab(Lab_m/2,Lcd_p+3);
    double Icd=Dcd(Lcd_m/2,Lab_p+3);
//    if(Lab_p>0 || Lcd_p>0 )
//    {
//        std::cout <<  "L+=" << Lab_p << ", " << Lcd_p << ", L-=" << Lab_m << ", " << Lcd_m << std::endl;
//        std::cout <<  afact << "*" << Iab << " + " << cfact << "*" << Icd << std::endl;
//    }
    return Pi12/8*(afact*Iab+cfact*Icd);
}

double SlaterIntegrals::DoExchangeSum(int la, int lb, int lc, int ld) const
{
    static Wigner3j w;
    assert(la==lc);
    assert(lb==ld);
    assert(la>=0);
    assert(lb>=0);
    int lmin=std::abs(la-lb);
    int lmax=la+lb;
    double ret=0.0;
    for (int l=lmin;l<=lmax;l+=2)
    {
        ret+=(*this)(l,la,lb,la,lb)*w(la,l,lb);
    }
    return 2*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
}


double SlaterIntegrals::D(double _a, int m,int n) const
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

double SlaterIntegrals::I(double _a, int m,int n) const
{
    assert(n>0);
    assert(m>=0);
    assert(m<=NMax);
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
        return Fact[m-1]*(n/(2*pow((a+c),m))+1/pow(_a,m));
}

double SlaterIntegrals::Dab(int m,int n) const
{
    return D(a,m,n);
}

double SlaterIntegrals::Dcd(int m,int n) const
{
    return D(c,m,n);
}
