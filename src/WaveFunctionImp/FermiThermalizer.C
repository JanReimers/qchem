// File: FermiThermalizer.C  Calculate level occupation factors at finit temperature.

#include "Imp/WaveFunction/FermiThermalizer.H"
#include "Imp/WaveFunction/ElectronDumper.H"
#include "oml/vector.h"
#include <iostream>
#include <cmath>
#include <cassert>

void fermi(const Vector<double>& e, const Vector<double>& g,
           double kT, double mu,
           double& f, double& fp);

double Ef(const Vector<double>& e, const Vector<double>& g, double& NumE);

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
    if (NumE>Sum(g))
    {
        std::cerr << "Too many electrons " << NumE << ", not enough levels " << Sum(g) << std::endl;
    }
    itsMu=Ef(e,g,itsNf);
    if (itskT>0)
    {
        double f,fp,dmu=0;
        do
        {
            fermi(e,g,itskT,itsMu,f,fp);
            dmu=-(f-NumE)/fp;
            itsMu+=dmu;
        }
        while (fabs(dmu)>1e-14);
    }
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

double Ef(const Vector<double>& e, const Vector<double>& g, double& Nf)
{
    Vector<double>::const_iterator be(e.begin());
    Vector<double>::const_iterator bg(g.begin());
    if (Nf==0) return *be-1.0;
    int n=0;
    double ret=0;
    double N=Nf;
    while (N >=0 )
    {
        assert(*bg>0);
        Nf=N/(*bg);
        N-=*bg;
        ret=*be;
        bg++;
        be++;
        n++;
    };
    if (N==0) ret=(ret+*be)/2.0; //If HOMO is full then goto midpoint of gap.
//  cout << "Ef=" << ret << " Nf=" << Nf << " nlevel=" << n << " N=" << N << std::endl;
    return ret;
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

