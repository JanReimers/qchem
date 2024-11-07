// File: SlaterCD.C  4 electron Charge distribution of Slater orbitals. 

#include "Imp/Integrals/SlaterCD.H"
#include "Imp/Integrals/PascalTriangle.H"
#include "Imp/Integrals/Factorials.H"
#include "oml/vector.h"

//
//  Ranges:  
//    0 <= k <= 2LMax  in steps of 2
//    3 <= Lab_p=la+lb+3+k <= 4LMax+3
//    1 <= Lcd_m=lc+ld+1-k <= 2LMax+1
//    3 <= Lab_m=la+lb+1-k <= 4LMax+3
//    1 <= Lcd_p=lc+ld+3+k <= 2LMax+1
//
//  Build up the derivative look up tables.
//
 SlaterCD::SlaterCD(double _eab, double _ecd, size_t _LMax)
 : eab(_eab), ecd(_ecd), LMax(_LMax), Iab(0,2*LMax+1,3,4*LMax+3), Icd(0,2*LMax+1,3,4*LMax+3)
 {
    assert(Iab.GetLimits()==Icd.GetLimits());
    Fill(Iab,0.0);
    Fill(Icd,0.0);
    Vector<double> f(0,2*LMax,0.0);
    const PascalTriangle& c1(PascalTriangle::thePascalTriangle); //Binomial coefficients.
    double eabcd=eab+ecd;
    for (size_t L2:Iab.cols())
    {
        double fL2=qchem::Fact[L2-1]; //(L2-1)!
        for (auto ik:f.indices()) f(ik)=fk(eab,eabcd,ik,L2);
        Iab(0,L2)=fL2/(eab*pow(eabcd,L2)); //This is what gets differentiated.
        for (size_t ik=1;ik<=2*LMax+1;ik++)
            for (size_t jk=0;jk<=ik-1;jk++)
                Iab(ik,L2)+=c1(ik-1,jk)*Iab(jk,L2)*f(ik-1-jk);  
            
        for (auto ik:f.indices()) f(ik)=fk(ecd,eabcd,ik,L2);
        Icd(0,L2)=fL2/(ecd*pow(eabcd,L2)); //This is what gets differentiated.
        for (size_t ik=1;ik<=2*LMax+1;ik++)
            for (size_t jk=0;jk<=ik-1;jk++)
                Icd(ik,L2)+=c1(ik-1,jk)*Icd(jk,L2)*f(ik-1-jk);  
    }
        
 }
 
 double SlaterCD::fk(double a, double ab, int k,int n)
{
    assert(n>0);
    assert(k>=0);
    assert(k<=qchem::NMax);
    return qchem::Fact[k]*(n/pow(ab,k+1)+1/pow(a,k+1));
}

Vector<double> SlaterCD::Coulomb_Rk(int la,int lc) const
{
    Vector<double> ret(la+lc+1,0.0);
    int i=1;
    for (int k=0;k<=2*std::min(la,lc);k+=2)
    {
        int Lab_p=2*la+3+k; // first term r_1^2
        int Lcd_m=2*lc+1-k; // first term r_2
        int Lab_m=2*la+1-k; // second term r_1
        int Lcd_p=2*lc+3+k; // second term r_2^2
        //cout << la << " " << lc << " " << k << " " << Lab_p << " " << Lcd_p << endl;
        ret(i++)=(2*k+1)*(Iab(Lab_m,Lcd_p)+Icd(Lcd_m,Lab_p));
    }
    return ret;
}

Vector<double> SlaterCD::ExchangeRk(int la,int lb) const
{
    int kmin=std::abs(la-lb);
    int kmax=la+lb;
    int N=(kmax-kmin)/2+1;
    Vector<double> ret(N,0.0);
    int i=1;
    for (int k=kmin;k<=kmax;k+=2)
    {
        int Lab_p=la+lb+3+k; // first term r_1^2
        int Lcd_m=la+lb+1-k; // first term r_2
        int Lab_m=la+lb+1-k; // second term r_1
        int Lcd_p=la+lb+3+k; 
        ret(i++)=(Iab(Lab_m,Lcd_p)+Icd(Lcd_m,Lab_p)); //(2*k+1)???
    }
    return ret;
}


