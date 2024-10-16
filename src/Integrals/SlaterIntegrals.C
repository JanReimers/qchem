#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Factorials.H"
#include "Misc/DFTDefines.H"
#include <cassert>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

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
//    assert(Lab_p%2==0);
//    assert(Lab_m%2==0);
//    assert(Lcd_p%2==0);
//    assert(Lcd_m%2==0);
    assert(Lab_m>=0);
    assert(Lcd_m>=0);

    assert(Lab_p+1<=qchem::NMax);
    assert(Lcd_p+1<=qchem::NMax);
//    cout << "Lab_m, Lcd_m, Lab_p, Lcd_p =" << Lab_m << " " << Lcd_m << " " << Lab_p << " " << Lcd_p << endl; 
    double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
    double cfact=qchem::Fact[Lab_p];
    double Iab=Dab(Lab_m,Lcd_p+1);
    double Icd=Dcd(Lcd_m,Lab_p+1);
//    if(Lab_p>0 || Lcd_p>0 )
//    {
//        std::cout <<  "L+=" << Lab_p << ", " << Lcd_p << ", L-=" << Lab_m << ", " << Lcd_m << std::endl;
//        std::cout <<  afact << "*" << Iab << " + " << cfact << "*" << Icd << std::endl;
//    }
//    cout << "afact*Iab+cfact*Icd=" << afact*Iab+cfact*Icd << endl;
    return afact*Iab+cfact*Icd;
}
double SlaterRadialIntegrals::DoExchangeSum(      int la, int lb, int lc, int ld) const
{
//    if (la==1 && lb==1 && lc==1 && ld==1)
//        cout << "DoExchangeSum (" << la << "," << lb << "," << lc << "," << ld << ")" << endl;
//                    
    static Wigner3j w; //Returns the **square** of the 3j symbol.
    assert(la==lc);
    assert(lb==ld);
    assert(la>=0);
    assert(lb>=0);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    double ret=0.0;
    for (int k=kmin;k<=kmax;k+=2)
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
    //cout << "D(a,k,n) a=" << _a << " k=" << k << " n=" << n << endl; 
    double ret=0.0;
    double I0=1/(_a*pow(eab+ecd,n));
    //cout << "I0=" << I0 << endl;
    double I1,I2,I3,I4,I5,I6,f0,f1,f2,f3,f4,f5;
   
    if (k==0) return I0;
    f0=f(_a,0,n);
    I1=I0*f0;
    if (k==1) return I1;
    f1=f(_a,1,n);
    I2=I1*f0 + I0*f1;
    if (k==2) return I2;

    f2=f(_a,2,n);
    I3=I2*f0 + 2*I1*f1 + I0*f2;
    if (k==3) return I3;
    
    f3=f(_a,3,n);
    I4=I3*f0 + 3*I2*f1 + 3*I1*f2 + I0*f3;
    if (k==4) return I4;
    
    f4=f(_a,4,n);
    I5=I4*f0 + 4*I3*f1 + 6*I2*f2 + 4*I1*f3 + I0*f4;
    if (k==5) return I5;
    
    f5=f(_a,5,n);
    I6=I5*f0 + 5*I4*f1 + 10*I3*f2 + 10*I2*f3 + 5*I2*f4 + I0*f5;
    if (k==6) return I6;
    
    
    cout << "k=" << k << " is too large.  Keep differentiating I(a) !!!" << endl;
    assert(false);
    return -1;
}


double SlaterRadialIntegrals::f(double _a, int k,int n) const
{
    assert(n>0);
    assert(k>=0);
    assert(k<=qchem::NMax);
    assert(_a==eab || _a==ecd);
    return qchem::Fact[k]*(n/pow(eab+ecd,k+1)+1/pow(_a,k+1));
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
