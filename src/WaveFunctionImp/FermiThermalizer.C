// File: FermiThermalizer.C  Calculate level occupation factors at finit temperature.

#include "Imp/WaveFunction/FermiThermalizer.H"
#include "Imp/WaveFunction/ElectronDumper.H"
#include "oml/vector.h"
#include "oml/imp/subvector.h"
#include <iostream>
#include <cmath>
#include <cassert>

void fermi(const Vector<double>& e, const Vector<double>& g,
           double kT, double mu,
           double& f, double& fp);

// return <Ef,Nf,last_index>        
std::tuple<double,double,int> Ef(const Vector<double>& e, const Vector<double>& g, double NumE);

FermiThermalizer::FermiThermalizer(const optr_vector1<EnergyLevel*>& el,
                                   double kT, double NumE)
    : itskT(kT)
    , itsMu(0)  //Acts as Ef when kT=0;
    , itsNf(NumE)  //fractional occupation at Ef when kT=0;
{
    int NLevel=el.size();
    Vector<double> e(NLevel),g(NLevel); //base 1 indexing.
    for (int i=0; i<NLevel; i++)
    {
        e(i+1)=el[i]->GetEnergy    ();
        g(i+1)=el[i]->GetDegeneracy();
    }
    //std::cout << "e,g=" << e << " " << g << std::endl;
    if (NumE>Sum(g))
    {
        std::cerr << "Too many electrons " << NumE << ", not enough levels " << Sum(g) << std::endl;
    }
    int last_index;
    std::tie(itsMu,itsNf,last_index)=Ef(e,g,NumE);
    if (itskT>0)
    {
        Vector<double> e1(last_index+1);
        Vector<double> g1(last_index+1);
        for (auto i:e1.indices())
        {
            e1(i)=e(i);
            g1(i)=g(i);
        }
        
        double f,fp,dmu=0;
        do
        {
            fermi(e1,g1,itskT,itsMu,f,fp);
            dmu=-(f-NumE)/fp;
            itsMu+=dmu;
        }
        while (fabs(dmu)>1e-14);
    }
    std::cout << "FermiThermalizer kT,Mu,Nf=" << itskT << " " << itsMu << " " << itsNf << std::endl;
}

double FermiThermalizer::GetOccupation(double e)
{
    double ret=0.0;
    if (itskT>0)
    {
        ret = 1/(1+exp((e-itsMu)/itskT));
        if (ret<1e-6) ret=0.0;
    }
    else
    {
        if (e<=itsMu) ret=1.0;
//        if (fabs(e-itsMu)<0.000001) ret=itsNf;
    }
    return ret;
}

std::tuple<double,double,int>  Ef(const Vector<double>& e, const Vector<double>& g, double NumE)
{
    Vector<double>::const_iterator be(e.begin());
    Vector<double>::const_iterator bg(g.begin());
    assert(NumE>0.0);
//    if (NumE==0) return *be-1.0;
    int n=0;
    double Ef=0,Nf=0.0;
    double N=NumE;
    while (N >=0 )
    {
        assert(*bg>0);
        Nf=N/(*bg);
        N-=*bg;
        Ef=*be;
        bg++;
        be++;
        n++;
//        std::cout << N << " " << Nf << " " << Ef << std::endl;
    };
    if (N==0) Ef=(Ef+*be)/2.0; //If HOMO is full then goto midpoint of gap.
//  cout << "Ef=" << ret << " Nf=" << Nf << " nlevel=" << n << " N=" << N << std::endl;
    return std::make_tuple(Ef,Nf,n);
}



void fermi(const Vector<double>& e, const Vector<double>& g,
           double kT, double mu,
           double& f, double& fp)
{
    double onekT=1.0/kT;
    Vector<double> ex=exp((e-mu)*onekT);
    Vector<double> ex1=1.0/(ex+1.0);
    f=Dot(g,ex1);
    ex1=DirectMultiply(ex1,ex1);
    ex1=DirectMultiply(ex1,ex);
    fp=onekT*Dot(g,ex1);
}

