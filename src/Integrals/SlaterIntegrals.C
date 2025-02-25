#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Factorials.H"
#include "Imp/Integrals/PascalTriangle.H"
#include "Imp/Misc/DFTDefines.H"
#include "oml/vector.h"
#include <cassert>
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

//double SlaterRadialIntegrals::FourPi2=4*4*Pi*Pi;
//
// SlaterRadialIntegrals::SlaterRadialIntegrals(double _eab, double _ecd)
//    : eab(_eab), ecd(_ecd)
//{
//};
//    
//
//double SlaterRadialIntegrals::R(int k,int la, int lb, int lc, int ld) const
//{
//    int Lab_p=la+lb+2+k; // first term r_1^2
//    int Lcd_m=lc+ld+1-k; // first term r_2
//    int Lab_m=la+lb+1-k; // second term r_1
//    int Lcd_p=lc+ld+2+k; // second term r_2^2
//    assert(Lab_m>=0);
//    assert(Lcd_m>=0);
//    assert(Lab_p+1<=qchem::NMax);
//    assert(Lcd_p+1<=qchem::NMax);
//    double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
//    double cfact=qchem::Fact[Lab_p];
//    double Iab=D(eab,Lab_m,Lcd_p+1);
//    double Icd=D(ecd,Lcd_m,Lab_p+1);
//    //cout << " SlaterRI::Rk " <<  Lab_p << " " << Lcd_m << " " << Lab_m << " " << Lcd_p << " " << Iab << " " << Icd << " " << eab << " " << ecd << endl;
//    return afact*Iab+cfact*Icd;
//}
//
//double SlaterRadialIntegrals::Coulomb(int lab, int lcd) const
//{
//    int Lab_p=lab+2; // first term r_1^2
//    int Lcd_m=lcd+1; // first term r_2
//    int Lab_m=lab+1; // second term r_1
//    int Lcd_p=lcd+2; // second term r_2^2
//    assert(Lab_m>=0);
//    assert(Lcd_m>=0);
//    assert(Lab_p+1<=qchem::NMax);
//    assert(Lcd_p+1<=qchem::NMax);
//    double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
//    double cfact=qchem::Fact[Lab_p];
//    double Iab=D(eab,Lab_m,Lcd_p+1);
//    double Icd=D(ecd,Lcd_m,Lab_p+1);
//
//    return afact*Iab+cfact*Icd;
//}
//
//double SlaterRadialIntegrals::Coulomb(int la, int lb, int lc, int ld) const
//{
//    assert(la==lb);
//    assert(lc==ld);
//    return FourPi2*R(0,la,lb,lc,ld);
////    return (2*la+1)*(2*lc+1)*R(0,la,lb,lc,ld);
//}
//
//double SlaterRadialIntegrals::DoExchangeSum(int la, int lb, int lc, int ld) const
//{
//    if (la==ld && lb==lc && la!=lb) return DoExchangeSum(la,lb,ld,lc);
//    assert(la==lc);
//    assert(lb==ld);
//    assert(la>=0);
//    assert(lb>=0);
//    int kmin=std::abs(la-lb);
//    int kmax=la+lb;
//    double ret=0.0;
//    for (int k=kmin;k<=kmax;k+=2)
//    {
//        assert((k+la+lb)%2==0);
//        //ret+=R(k,la,lb,la,lb)*Wigner3j::theW3j(la,k,lb); //What about *(2k+1) ??
//        ret+=R(k,la,lb,la,lb)*AngularIntegrals::Exchange(k,la,lb); //What about *(2k+1) ??
//
//    }
//    return FourPi2*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
////    return ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
//}
//
//double SlaterRadialIntegrals::Coulomb(int la, int lb, int lc, int ld,int ma, int mb, int mc, int md) const
//{
//    assert(la==lb);
//    assert(lc==ld);
//    assert(la>=0);
//    assert(lc>=0);
//    assert(ma>=-la);
//    assert(ma<= la);
//    assert(mc>=-lc);
//    assert(mc<= lc);
//    int kmin=0;
//    int kmax=2*std::min(la,lc);
//    double ret=0.0;
//    for (int k=kmin;k<=kmax;k+=2)
//    {
//        //cout << "Rk*(2k+1) , A = " << R(k,la,la,lc,lc)*(2*k+1) << " " << AngularIntegrals::Coulomb(k,la,lc,ma,mc) << endl;
//        ret+=R(k,la,la,lc,lc)*(2*k+1)*AngularIntegrals::Coulomb(k,la,lc,ma,mc); //What about *(2k+1) ??
//    }
//    return FourPi2*(2*la+1)*(2*lc+1)*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
////    return ret; //Compensate for
//}
//
//double SlaterRadialIntegrals::DoExchangeSum(int la, int lb, int lc, int ld, int ma, int mb, int mc, int md) const
//{
//    if (la==ld && lb==lc && la!=lb) 
//        return DoExchangeSum(la,lb,ld,lc,ma,mb,md,mc);
//   
//    assert(la==lc);
//    assert(lb==ld);
////    if (la==lb && ma==md && mb==mc && ma!=mc && mb!=md) 
////        return DoExchangeSum(la,lb,ld,lc,ma,mb,md,mc);
////    assert(ma==mc);
////    assert(mb==md);
//    assert(la>=0);
//    assert(lb>=0);
//    assert(ma>=-la);
//    assert(ma<= la);
//    assert(mb>=-lb);
//    assert(mb<= lb);
//    int kmin=std::abs(la-lb);
//    int kmax=la+lb;
//    double ret=0.0;
//    for (int k=kmin;k<=kmax;k+=2)
//    {
//        assert((k+la+lb)%2==0);
//       //cout << "Rk*(2k+1) , A = " << R(k,la,la,lc,lc)*(2*k+1) << " " << AngularIntegrals::Coulomb(k,la,lc,ma,mc) << endl;
//
//        ret+=R(k,la,lb,la,lb)*AngularIntegrals::Exchange(k,la,lb,ma,mb); //What about *(2k+1) ??
//    }
//    return FourPi2*(2*la+1)*(2*lb+1)*ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
////    return ret; //Compensate for factor if 1/2 built into the Wigner3j lookup tables.
//}
//
////
////  Auto differentiate I_0 k times.
////
//double SlaterRadialIntegrals::D(double _a, int k,int n) const
//{
//    Vector<double> I(0,k,0.0),f(0,k-1,0.0);
//    const PascalTriangle& c1(PascalTriangle::thePascalTriangle);
//    I(0)=1/(_a*pow(eab+ecd,n));
//    for (auto ik:f.indices()) f(ik)=fk(_a,ik,n);
//    for (int ik=1;ik<=k;ik++)
//         for (int jk=0;jk<=ik-1;jk++)
//            I(ik)+=c1(ik-1,jk)*I(jk)*f(ik-1-jk);             
//    return I(k);
//}
//
//double SlaterRadialIntegrals::fk(double _a, int k,int n) const
//{
//    assert(n>0);
//    assert(k>=0);
//    assert(k<=qchem::NMax);
//    assert(_a==eab || _a==ecd);
//    return qchem::Fact[k]*(n/pow(eab+ecd,k+1)+1/pow(_a,k+1));
//}
////
////
////
//Vector<double> Rk(int la,int lc,double eab,double ecd)
//{
//    Vector<double> ret(la+lc+1);
//    int i=1;
//    int kmax=2*std::min(la,lc);
//    for (int k=0;k<=kmax;k+=2)
//    {
//        int Lab_p=2*la+2+k; // first term r_1^2
//        int Lcd_m=2*lc+1-k; // first term r_2
//        int Lab_m=2+la+1-k; // second term r_1
//        int Lcd_p=2*lc+2+k; // second term r_2^2
//        assert(Lab_m>=0);
//        assert(Lcd_m>=0);
//        assert(Lab_p+1<=qchem::NMax);
//        assert(Lcd_p+1<=qchem::NMax);
//        double afact=qchem::Fact[Lcd_p]; //These ab and cd are reversed on purpose.
//        double cfact=qchem::Fact[Lab_p];
//        double Iab=D(eab,ecd,Lab_m,Lcd_p+1);
//        double Icd=D(ecd,eab,Lcd_m,Lab_p+1);
//        ret(i++)=afact*Iab+cfact*Icd;
//    }
//    return ret;
//}
