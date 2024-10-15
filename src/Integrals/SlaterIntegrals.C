#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Factorials.H"
#include "Misc/DFTDefines.H"
#include <cassert>
#include <cmath>
#include <iostream>


 SlaterRadialIntegrals::SlaterRadialIntegrals(double _eab, double _ecd)
    : eab(_eab), ecd(_ecd)
{
    if (qchem::DFact[0]!=1.0) qchem::InitFactorials();
};
    

double SlaterRadialIntegrals::operator()   (int k,int la, int lb, int lc, int ld) const
{
    int Lab_p=la+lb+2+k; // first term r_1^2
    int Lcd_m=lc+ld+1-k; // first term r_2
    int Lab_m=la+lb+1-k; // second term r_1
    int Lcd_p=lc+ld+2+k; // second term r_2^2
    //Check that everything is even.
    assert(Lab_p%2==0);
    assert(Lab_m%2==0);
    assert(Lcd_p%2==0);
    assert(Lcd_m%2==0);
    assert(Lab_m>=0);
    assert(Lcd_m>=0);

    assert(Lab_p+1<=qchem::NMax);
    assert(Lcd_p+1<=qchem::NMax);
    double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
    double cfact=qchem::Fact[Lab_p];
    double Iab=Dab(Lab_m,Lcd_p+1);
    double Icd=Dcd(Lcd_m,Lab_p+1);
//    if(Lab_p>0 || Lcd_p>0 )
//    {
//        std::cout <<  "L+=" << Lab_p << ", " << Lcd_p << ", L-=" << Lab_m << ", " << Lcd_m << std::endl;
//        std::cout <<  afact << "*" << Iab << " + " << cfact << "*" << Icd << std::endl;
//    }
    return afact*Iab+cfact*Icd;
}
double SlaterRadialIntegrals::DoExchangeSum(      int la, int lb, int lc, int ld) const
{
    static Wigner3j w; //Returns the **square** of the 3j symbol.
    assert(la==lc);
    assert(lb==ld);
    assert(la>=0);
    assert(lb>=0);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    double ret=0.0;
    for (int k=kmin;k<=kmax;k++)
    {
        ret+=(*this)(k,la,lb,la,lb)*w(la,k,lb);
    }
    return 2*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
}

double SlaterRadialIntegrals::Dab(int m,int n) const
{
    return D(eab,m,n);
}

double SlaterRadialIntegrals::Dcd(int m,int n) const
{
    return D(ecd,m,n);
}

// _a can be alpha or beta from the tech. notes.
double SlaterRadialIntegrals::D(double _a, int k,int n) const
{
    double ret=0.0;
    switch (k)
    {
    case 0:
    {
        ret=I(_a,0,n);
        break;
    }
    case 1:
    {
        ret=I(_a,0,n)*I(_a,1,n); //I * I1
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


double SlaterRadialIntegrals::I(double _a, int k,int n) const
{
    assert(n>0);
    assert(k>=0);
    assert(k<=qchem::NMax);
    assert(_a==eab || _a==ecd);
    if (k==0)
        return 1.0/(_a*pow(eab+ecd,n));
    else
        return n/pow(eab+ecd,k)+1/pow(_a,k);
}


//
//
//double SlaterRadialIntegrals::R(int k) const
//{
//    return qchem::Fact[nac-k-1]/intpow(eac,nac-k)*Habcd(k)
//         + qchem::Fact[nbd-k-1]/intpow(ebd,nbd-k)*Hbadc(k);
//}
//
//double SlaterRadialIntegrals::Habcd(int k) const
//{
//   double  Hk=0.0;
//   for (int i=0;i<nac-k-1;i++) Hk+=qchem::Fact[nbd+k+i]*intpow(eac,i)/intpow(eabcd,nbd+k+i+1)/qchem::Fact[i];
//   return Hk;
//}
//double SlaterRadialIntegrals::Hbadc(int k) const
//{
//   double  Hk=0.0;
//   for (int i=0;i<nac-k-1;i++) Hk+=qchem::Fact[nac+k+i]*intpow(ebd,i)/intpow(eabcd,nac+k+i+1)/qchem::Fact[i];
//   return Hk;
//}
//
