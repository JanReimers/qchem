#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Factorials.H"
#include "Misc/DFTDefines.H"
#include "oml/vector.h"
#include <cassert>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

 SlaterRadialIntegrals::SlaterRadialIntegrals(double _eab, double _ecd)
    : eab(_eab), ecd(_ecd)
{
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
//    double Iab=Dab(Lab_m,Lcd_p+1);
//    double Icd=Dcd(Lcd_m,Lab_p+1);
    
    double Iab=D1(eab,Lab_m,Lcd_p+1);
    double Icd=D1(ecd,Lcd_m,Lab_p+1);
    
//    double errab=fabs(Iab-Iab1)/Iab;
//    if (errab>=1e-13)
//    {
//        cout << Lab_m << " " << Lcd_p+1 << endl;
//        cout << Iab << " " << Iab1 << " " << errab << endl;
//    }
//    assert(errab<1e-13);

//    if(Lab_p>0 || Lcd_p>0 )
//    {
//        std::cout <<  "L+=" << Lab_p << ", " << Lcd_p << ", L-=" << Lab_m << ", " << Lcd_m << std::endl;
//        std::cout <<  afact << "*" << Iab << " + " << cfact << "*" << Icd << std::endl;
//    }
//    cout << "afact*Iab+cfact*Icd=" << afact*Iab+cfact*Icd << endl;
    return afact*Iab+cfact*Icd;
}

double SlaterRadialIntegrals::operator()   (int lab, int lcd) const
{
    int Lab_p=lab+2; // first term r_1^2
    int Lcd_m=lcd+1; // first term r_2
    int Lab_m=lab+1; // second term r_1
    int Lcd_p=lcd+2; // second term r_2^2
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
    assert(la==lc);
    assert(lb==ld);
    assert(la>=0);
    assert(lb>=0);
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    double ret=0.0;
    for (int k=kmin;k<=kmax;k+=2)
    {
        ret+=(*this)(k,la,lb,la,lb)*Wigner3j::theW3j(la,k,lb);
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
    double I0=1/(_a*pow(eab+ecd,n));
    //cout << "I0=" << I0 << endl;
    double I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
   
    if (k==0) return I0;
    f0=fk(_a,0,n);
    I1=I0*f0;
    if (k==1) return I1;
    f1=fk(_a,1,n);
    I2=I1*f0 + I0*f1;
    if (k==2) return I2;

    f2=fk(_a,2,n);
    I3=I2*f0 + 2*I1*f1 + I0*f2;
    if (k==3) return I3;
    
    f3=fk(_a,3,n);
    I4=I3*f0 + 3*I2*f1 + 3*I1*f2 + I0*f3;
    if (k==4) return I4;
    
    f4=fk(_a,4,n);
    I5=I4*f0 + 4*I3*f1 + 6*I2*f2 + 4*I1*f3 + I0*f4;
    if (k==5) return I5;
    
    f5=fk(_a,5,n);
    I6=I5*f0 + 5*I4*f1 + 10*I3*f2 + 10*I2*f3 + 5*I1*f4 + I0*f5;
    if (k==6) return I6;
    
    f6=fk(_a,6,n);
    I7=I6*f0 + 6*I5*f1 + 15*I4*f2 + 20*I3*f3 + 15*I2*f4 + 6*I1*f5 + I0*f6;
    if (k==7) return I7;

    f7=fk(_a,7,n);
    I8=I7*f0 + 7*I6*f1 + 21*I5*f2 + 35*I4*f3 + 35*I3*f4 + 21*I2*f5 + 7*I1*f6 + I0*f7;
    if (k==8) return I8;
    
    f8=fk(_a,8,n);
    I9=I8*f0 + 8*I7*f1 + 28*I6*f2 + 56*I5*f3 + 70*I4*f4 + 56*I3*f5 + 28*I2*f6 + 8*I1*f7 + I0*f8;
    if (k==9) return I9;

    f9=fk(_a,9,n);
    I10=I9*f0 + 9*I8*f1 + 36*I7*f2 + 82*I6*f3 + 126*I5*f4 + 126*I4*f5 + 82*I3*f6 + 36*I2*f7 + 9*I1*f8+ I0*f9;
    if (k==10) return I10;

    cout << "k=" << k << " is too large.  Keep differentiating I(a) !!!" << endl;
    assert(false);
    return -1;
}

double SlaterRadialIntegrals::D1(double _a, int k,int n) const
{
    Vector<double> c(0,k-1,1.0),I(0,k,0.0),f(0,k-1,0.0);
//    Fill(c,1.0);
//    Fill(I,0.0);
    I(0)=1/(_a*pow(eab+ecd,n));
    for (auto ik:f.indices()) f(ik)=fk(_a,ik,n);
    for (int ik=1;ik<=k;ik++)
    {
        // build next row in Pascal's triangle.
         double ctmp=1.0,ctmp1;
         for (int jk=1;jk<ik-1;jk++)
         {
            ctmp1=c(jk);
            c(jk)=ctmp+c(jk);
            ctmp=ctmp1;
         }
         for (int jk=0;jk<=ik-1;jk++)
            I(ik)+=c(jk)*I(jk)*f(ik-1-jk);             
    }
    return I(k);
}

double SlaterRadialIntegrals::fk(double _a, int k,int n) const
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
