// File: Gaussian::Rk.C  2 electron Charge distribution of Gaussian orbitals. 
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
module qchem.BasisSet.Atom.radial.GaussianRk;
import qchem.BasisSet.Atom.radial.PascalTriangle;
import Common.Constants;
import Common.Factorials;
import oml;

using std::cout;
using std::endl;

namespace Gaussian
{
//
//  Ranges:  
//    0 <= k <= 2LMax  in steps of 2
//    1 <= Lab_p=0.5*(la+lb+2+k) <= 2LMax+1
//    0 <= Lcd_m=0.5*(lc+ld  -k) <=  LMax
//    0 <= Lab_m=0.5*(la+lb  -k) <= 2LMax+1
//    1 <= Lcd_p=0.5*(lc+ld+2+k) <=  LMax
//
//  Build up the derivative look up tables.
//
RkEngine::RkEngine(double _eab, double _ecd, size_t _LMax)
 : eab(_eab), ecd(_ecd), LMax(_LMax), Iab(0,LMax,1,2*LMax+1), Icd(0,LMax,1,2*LMax+1)
 {
 //   cout << "RkEngine eab,ecd,LMax=" << eab << " " << ecd << " " << LMax << endl;
    assert(Iab.GetLimits()==Icd.GetLimits());
    Fill(Iab,0.0);
    Fill(Icd,0.0);
    Vector<double> f(0,LMax,0.0);
    const PascalTriangle& c1(PascalTriangle::thePascalTriangle); //Binomial coefficients.
    double eabcd=eab+ecd;
    
    for (size_t L2:Iab.cols())
    {
        double fL2=qchem::DFact[2*L2-1]/pow(2,L2-1); //sqrt(pi)*(2*n-1)!!/2^n/4
        for (auto ik:f.indices()) f(ik)=fk(eab,eabcd,ik,L2);
        Iab(0,L2)=fL2/(eab*pow(eabcd,L2+0.5)); //This is what gets differentiated.
        //cout << "L2,Iab(0,L2) " << L2 << " " << Iab(0,L2) << endl;
        for (size_t ik=1;ik<=LMax;ik++)
            for (size_t jk=0;jk<=ik-1;jk++)
                Iab(ik,L2)+=c1(ik-1,jk)*Iab(jk,L2)*f(ik-1-jk);  
            
        for (auto ik:f.indices()) f(ik)=fk(ecd,eabcd,ik,L2);
        Icd(0,L2)=fL2/(ecd*pow(eabcd,L2+0.5)); //This is what gets differentiated.
        //cout << "L2,Icd(0,L2) " << L2 << " " << Icd(0,L2) << endl;
        for (size_t ik=1;ik<=LMax;ik++)
            for (size_t jk=0;jk<=ik-1;jk++)
            {
                Icd(ik,L2)+=c1(ik-1,jk)*Icd(jk,L2)*f(ik-1-jk);  
//                if (ik==1 && L2==1)
//                {
//                    cout << jk << " " << c1(ik-1,jk) << " " << Icd(jk,L2) << " " << f(ik-1-jk) << endl;
//                }
                
            }
    }
        
 }
 
 double RkEngine::fk(double a, double ab, size_t k,size_t n)
{
    assert(n>0);
    assert(k>=0);
    assert(k<=qchem::NMax);
    return qchem::Fact[k]*((n+0.5)/pow(ab,k+1)+1/pow(a,k+1));
}

double RkEngine::Coulomb_R0(size_t la,size_t lc) const
{
    assert(la<=LMax);
    assert(lc<=LMax);
    size_t Lab_p=la+1; // (la+lb+2)/2 
    size_t Lcd_m=lc;   // (lc+ld)/2   
    size_t Lab_m=la;   // (la+lb)/2   
    size_t Lcd_p=lc+1; // (lc+ld+2)/2
    //cout << "Lab_m Lcd_m Lab_p Lcd_p" << Lab_m << " " << Lcd_m << " " << Lab_p << " " << Lcd_p << endl;
    //cout << "Iab Icd = " << Iab(Lab_m,Lcd_p) << " " << Icd(Lcd_m,Lab_p) << endl;
    return Pi12/8*(Iab(Lab_m,Lcd_p)+Icd(Lcd_m,Lab_p));
}

Vector<double> RkEngine::Coulomb_Rk(size_t la,size_t lc) const
{
    assert(la<=LMax);
    assert(lc<=LMax);
    Vector<double> ret(la+lc+1,0.0);
    size_t i=1;
    for (size_t k=0;k<=2*std::min(la,lc);k+=2)
    {
        size_t Lab_p=la+1+k/2; 
        size_t Lcd_m=lc-k/2; 
        size_t Lab_m=la-k/2; 
        size_t Lcd_p=lc+1+k/2;
        //cout << la << " " << lc << " " << k << " " << Lab_p << " " << Lcd_p << endl;
        ret(i++)=Pi12/8*(Iab(Lab_m,Lcd_p)+Icd(Lcd_m,Lab_p));
    }
    return ret;
}

Vector<double> RkEngine::ExchangeRk(size_t la,size_t lb) const
{
    assert(la<=LMax);
    assert(lb<=LMax);
    size_t kmin=std::abs((int)la-(int)lb);
    size_t kmax=la+lb;
    size_t N=(kmax-kmin)/2+1;
    Vector<double> ret(N,0.0);
    size_t i=1;
    for (size_t k=kmin;k<=kmax;k+=2)
    {
        assert((la+lb+k)%2==0);
        size_t Lab_p=(la+lb+k)/2+1; 
        size_t Lcd_m=(la+lb-k)/2; 
        size_t Lab_m=(la+lb-k)/2; 
        size_t Lcd_p=(la+lb+k)/2+1;

        ret(i++)=Pi12/8*(Iab(Lab_m,Lcd_p)+Icd(Lcd_m,Lab_p)); //(2*k+1)???
    }
    return ret;
}

} //namespace

