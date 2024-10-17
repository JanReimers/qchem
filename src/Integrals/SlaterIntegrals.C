#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Factorials.H"
#include "Imp/Integrals/PascalTriangle.H"
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
    assert(Lab_m>=0);
    assert(Lcd_m>=0);
    assert(Lab_p+1<=qchem::NMax);
    assert(Lcd_p+1<=qchem::NMax);
    double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
    double cfact=qchem::Fact[Lab_p];
    double Iab=D(eab,Lab_m,Lcd_p+1);
    double Icd=D(ecd,Lcd_m,Lab_p+1);
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
    double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
    double cfact=qchem::Fact[Lab_p];
    double Iab=D(eab,Lab_m,Lcd_p+1);
    double Icd=D(ecd,Lcd_m,Lab_p+1);
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
        ret+=(*this)(k,la,lb,la,lb)*Wigner3j::theW3j(la,k,lb);
    return 2*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
}

//
//  Auto differentiate I_0 k times.
//
double SlaterRadialIntegrals::D(double _a, int k,int n) const
{
    Vector<double> I(0,k,0.0),f(0,k-1,0.0);
    const PascalTriangle& c1(PascalTriangle::thePascalTriangle);
    I(0)=1/(_a*pow(eab+ecd,n));
    for (auto ik:f.indices()) f(ik)=fk(_a,ik,n);
    for (int ik=1;ik<=k;ik++)
         for (int jk=0;jk<=ik-1;jk++)
            I(ik)+=c1(ik-1,jk)*I(jk)*f(ik-1-jk);             
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

