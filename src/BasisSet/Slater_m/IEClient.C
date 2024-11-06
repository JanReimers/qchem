
#include "Imp/BasisSet/Slater_m/IEClient.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/PascalTriangle.H"

using std::cout;
using std::endl;

template <class T> void FillPower(Vector<T>& arr,T start, T stop);
//template void FillPower(Vector<double>& arr,double start, double stop);

namespace Slater_m
{
    
void IrrepIEClient::Init(double minexp,double maxexp,size_t L, int m)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      Fill(Ms,m);
      for (auto i:es.indices())  ns(i)=SlaterNorm(es(i),Ns(i));
}

void IEClient::Append(const IrrepIEClient* ic)
{
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ns.SetLimits(N,true);
    Ls.SetLimits(N,true);
    Ms.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ns(j)=ic->Ns(i);
        Ls(j)=ic->Ls(i);
        Ms(j)=ic->Ms(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
        
        size_t index=unique_es.size();
        if (const auto &ie =unique_es.find(es(j));ie==unique_es.end())
            unique_es[es(j)]=index;
        else 
            index=ie->second;
        
        es_indices.push_back(index);
        
        if (const auto &il =L_indices.find(Ls(j));il==L_indices.end())
            L_indices[Ls(j)]=std::vector<size_t>();
        
        L_indices[Ls(j)].push_back(j);
        
        
    }
//    for (auto e:unique_es) cout << e.first << " ";
//    cout << endl;
//    for (auto i:es_indexes) cout << i << " ";
//    cout << endl;

}

const std::vector<size_t>& IEClient::indices(size_t l) const
{
    auto i=L_indices.find(l);
    assert(i!=L_indices.end());
    return i->second;
}

using std::cout;
using std::endl;


//
//  Ranges: All 
//    0 <= k <= 2LMax  in steps of 2
//    3 <= Lab_p=la+lb+3+k <= 2*(LMax+1) + 2LMax
//    1 <= Lcd_m=lc+ld+1-k <= 2LMax+1
//    3 <= Lab_m=la+lb+1-k <= 2*(LMax+1) + 2LMax
//    1 <= Lcd_p=lc+ld+3+k <= 2LMax+1
//
 SlaterCD::SlaterCD(double _eab, double _ecd, size_t _LMax)
 : eab(_eab), ecd(_ecd), LMax(_LMax), Iab(0,2*LMax+1,3,4*LMax+3), Icd(0,2*LMax+1,3,4*LMax+3)
 {
    assert(Iab.GetLimits()==Icd.GetLimits());
    Fill(Iab,0.0);
    Fill(Icd,0.0);
    Vector<double> f(0,2*LMax,0.0);
    const PascalTriangle& c1(PascalTriangle::thePascalTriangle);
    double eabcd=eab+ecd;
    for (size_t L2:Iab.cols())
    {
        double fL2=qchem::Fact[L2-1]; //(L2-1)!
        for (auto ik:f.indices()) f(ik)=fk(eab,eabcd,ik,L2);
        Iab(0,L2)=fL2/(eab*pow(eabcd,L2));
        for (size_t ik=1;ik<=2*LMax+1;ik++)
            for (size_t jk=0;jk<=ik-1;jk++)
                Iab(ik,L2)+=c1(ik-1,jk)*Iab(jk,L2)*f(ik-1-jk);  
            
        for (auto ik:f.indices()) f(ik)=fk(ecd,eabcd,ik,L2);
        Icd(0,L2)=fL2/(ecd*pow(eabcd,L2));
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

//double SlaterCD::D(double a, double ab, int k,int n) 
//{
//    Vector<double> I(0,k,0.0),f(0,k-1,0.0);
//    const PascalTriangle& c1(PascalTriangle::thePascalTriangle);
//    I(0)=1/(a*pow(ab,n));
//    for (auto ik:f.indices()) f(ik)=fk(a,ab,ik,n);
//    for (int ik=1;ik<=k;ik++)
//         for (int jk=0;jk<=ik-1;jk++)
//            I(ik)+=c1(ik-1,jk)*I(jk)*f(ik-1-jk);             
//    return I(k);
//}

//double SlaterCD::R(int k,int la, int lb, int lc, int ld) const
//{
//    int Lab_p=la+lb+3+k; // first term r_1^2
//    int Lcd_m=lc+ld+1-k; // first term r_2
//    int Lab_m=la+lb+1-k; // second term r_1
//    int Lcd_p=lc+ld+3+k; // second term r_2^2
//    assert(Lab_m>=0);
//    assert(Lcd_m>=0);
//    assert(Lab_p+1<=qchem::NMax);
//    assert(Lcd_p+1<=qchem::NMax);
//    double afact=qchem::Fact[Lcd_p-1]; //These ab and cd are reversed on purpose.
//    double cfact=qchem::Fact[Lab_p-1];
//    double Iab=D(eab,eab+ecd,Lab_m,Lcd_p);
//    double Icd=D(ecd,eab+ecd,Lcd_m,Lab_p);
//    return afact*Iab+cfact*Icd;
//}

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



} //namespace
