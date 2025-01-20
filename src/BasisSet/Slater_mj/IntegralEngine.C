// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater_mj/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_mj/IEClient.H" 
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"
#include <iomanip>
#include <Imp/Integrals/SlaterCD.H>
#include <iostream>
#include <Imp/Symmetry/OkmjQN.H>

using std::cout;
using std::endl;

namespace Slater_mj
{
    
DiracIntegralEngine::DiracIntegralEngine()
    : itsLargeIE(new Slater_m::IntegralEngine())
    , itsSmallIE(new Small_IntegralEngine())
    {
     assert(itsLargeIE);
     assert(itsSmallIE);
    }

const Dirac_IrrepIEClient* DiracIntegralEngine::dcast(iec_t* iec)
{
    const Dirac_IrrepIEClient* diec=dynamic_cast<const Dirac_IrrepIEClient*>(iec);
    assert(diec);
    return diec;
}

DiracIntegralEngine::SMat DiracIntegralEngine::merge_diag(const SMat& l,const SMat& s)
{
    size_t Nl=l.GetNumRows();
    size_t Ns=s.GetNumRows();
    SMat ls(Nl+Ns);
    Fill(ls,0.0);
    for (auto i:l.rows())
        for (auto j:l.cols(i))
            ls(i,j)=l(i,j);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            ls(Nl+i,Nl+j)=s(i,j);
    return ls;
}

DiracIntegralEngine::SMat DiracIntegralEngine::merge_off_diag(const Mat& ls)
{
    size_t Nl=ls.GetNumRows();
    size_t Ns=ls.GetNumCols();
    assert(Nl==Ns);
    SMat k(Nl+Ns);
    Fill(k,0.0);
    for (auto i:ls.rows())
        for (auto j:ls.cols())
            k(i,Ns+j)=ls(i,j);
   
    return k;
}

DiracIntegralEngine::RVec DiracIntegralEngine::merge(const RVec& l,const RVec& s)
{
    size_t Nl=l.size();
    size_t Ns=s.size();
    RVec ls(Nl+Ns);
    Fill(ls,0.0);
    int i=1;
    for (auto c:l) ls(i++)=c;
    for (auto c:s) ls(i++)=c;
    return ls;
}

ERI4 DiracIntegralEngine::merge(const ERI4& LL,const ERI4& LS,const ERI4& SL,const ERI4& SS)
{
    size_t Nl=LL.size();
    size_t Ns=SS.size();
    ERI4 J(Nl+Ns,Nl+Ns);
    for (auto i:LL.rows())
        for (auto j:LL.cols(i))
            J(i,j)=LL(i,j);
    for (auto i:LS.rows())
        for (auto j:LS.cols(i))
            J(i,Nl+j)=LS(i,j);
    for (auto i:SL.rows())
        for (auto j:SL.cols(i))
            J(Nl+i,j)=SL(i,j);
    for (auto i:SS.rows())
        for (auto j:SS.cols(i))
            J(Nl+i,Nl+j)=SS(i,j);
    return J;
}

DiracIntegralEngine::SMat DiracIntegralEngine::MakeOverlap  (iec_t* a) const
{
    auto da=dcast(a);
    SMat ol=itsLargeIE->MakeOverlap(da->itsLargeIEC);
    SMat os=itsSmallIE->MakeOverlap(da->itsSmallIEC);
    return merge_diag(ol,os);
}

DiracIntegralEngine::SMat DiracIntegralEngine::MakeKinetic  (iec_t* a) const
{
    auto da=dcast(a);
    Mat kls=-2.0*itsLargeIE->MakeKinetic(da->itsLargeIEC,da->itsSmallIEC);
    return merge_off_diag(kls);
}

DiracIntegralEngine::Mat DiracIntegralEngine::MakeKinetic(iec_t* a,iec_t* b) const
{
    assert(false);
    return Mat();
}

DiracIntegralEngine::SMat DiracIntegralEngine::MakeNuclear  (iec_t* a, const Cluster& cl) const
{
    auto da=dcast(a);
    SMat vl=itsLargeIE->MakeNuclear(da->itsLargeIEC,cl);
    SMat vs=itsSmallIE->MakeNuclear(da->itsSmallIEC,cl);
    return merge_diag(vl,vs);
}
DiracIntegralEngine::SMat DiracIntegralEngine::MakeRestMass(iec_t* a) const
{
    static const double f=-2.0*c_light*c_light;
    auto da=dcast(a);
    SMat rl(da->itsLargeIEC->size());
    Fill(rl,0.0);
    SMat rs=f*itsSmallIE->MakeOverlap(da->itsSmallIEC);
    return merge_diag(rl,rs);
}

DiracIntegralEngine::RVec DiracIntegralEngine::MakeCharge  (iec_t* a) const
{
    auto da=dcast(a);
    RVec cl=itsLargeIE->MakeCharge(da->itsLargeIEC);
    RVec cs=itsSmallIE->MakeCharge(da->itsSmallIEC);
    return merge(cl,cs);
}

DiracIntegralEngine::SMat DiracIntegralEngine::MakeRepulsion  (iec_t* a) const
{
    assert(false);
    return SMat();
}
DiracIntegralEngine::Mat DiracIntegralEngine::MakeRepulsion  (iec_t*,iec_t*) const
{
    assert(false);
    return Mat();
}
DiracIntegralEngine::ERI3 DiracIntegralEngine::MakeOverlap3C  (iec_t* ab,iec_t* c) const
{
    assert(false);
    return ERI3();
}
DiracIntegralEngine::ERI3 DiracIntegralEngine::MakeRepulsion3C(iec_t* ab,iec_t* c) const
{
    assert(false);
    return ERI3();
}

ERI4 DiracIntegralEngine::MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const
{
    auto da=dcast(a);
    auto dc=dcast(c);
    assert(da->itsLargeIEC->kappa==-da->itsSmallIEC->kappa);
    assert(dc->itsLargeIEC->kappa==-dc->itsSmallIEC->kappa);
    ERI4 JLL=itsLargeIE->MakeDirect(da->itsLargeIEC,dc->itsLargeIEC);
    ERI4 JLS=itsSmallIE->MakeDirect(da->itsLargeIEC,dc->itsSmallIEC);
    ERI4 JSL=itsSmallIE->MakeDirect(da->itsSmallIEC,dc->itsLargeIEC);
    ERI4 JSS=itsSmallIE->MakeDirect(da->itsSmallIEC,dc->itsSmallIEC);
    
    return merge(JLL,JLS,JSL,JSS); 
}
ERI4 DiracIntegralEngine::MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const
{
    assert(false);
    return ERI4();    
}
DiracIntegralEngine::RVec DiracIntegralEngine::Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const
{
    assert(false);
    return RVec();    
}
DiracIntegralEngine::RVec DiracIntegralEngine::ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const
{
    assert(false);
    return RVec();    
}
Vector<double> DiracIntegralEngine::loop_4_direct  (size_t id, size_t la, size_t lc) const
{
   assert(false);
    return Vector<double>();        
}
Vector<double> DiracIntegralEngine::loop_4_exchange(size_t id, size_t la, size_t lc) const
{
    assert(false);
    return Vector<double>();    
}

const IrrepIEClient* Small_IntegralEngine::dcast(iec_t* iec)
{
    const IrrepIEClient* diec=dynamic_cast<const IrrepIEClient*>(iec);
    assert(diec);
    return diec;
}
// Overlap gets called with l*2
double Small_IntegralEngine::Overlap  (double ea , double eb,size_t l2) const
{
    assert(l2%2==0);
    return 2.0*Slater_m::IntegralEngine::Kinetic(ea,eb,l2/2); //Kinetic already has 4*Pi
}

//  This is anew one <a|p^2/r|b> !
//  For this one we actually need to know the sign kappa = l or -l-1
double Small_IntegralEngine::Nuclear(double ea, double eb,size_t l) const
{
    assert(l==0);
    int kappa = -l -1;
    return ea*eb*SlaterIntegral(ea+eb,-2*kappa-1);
   
}

//
//  Three cases to deal with:
//  1) a=L, c=S
//  2) a=S, c=L
//  3) a=S, c=S
//
ERI4 Small_IntegralEngine::MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const
{
    auto da=dcast(a);
    auto dc=dcast(c);  
    if (da->Large() && !dc->Large())
        return MakeDirectLS(da,dc);
    else if (!da->Large() && dc->Large())
        return MakeDirectSL(da,dc);
    else if (!da->Large() && !dc->Large())
        return MakeDirectSS(da,dc);
    assert(false);
    return ERI4();
}

ERI4 Small_IntegralEngine::MakeDirectLS(const IrrepIEClient* a, const IrrepIEClient* c) const
{
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        double ea=a->es_indices[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            double ec=c->es_indices[ic-1];
            //loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            // Small sector kappa is negative of the large sector kappa, which shifts lc but obly for the angular integrals.
            RVec Akac=Coulomb_AngularIntegrals(la,Omega_kQN::l(-c->kappa),0,0);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                double eb=a->es_indices[ib-1];
                //loop_3(a->es_indices[ib-1]);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << Jab(ic,id) << endl;    
                        assert(false);
                    }
                    double ed=c->es_indices[id-1];
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    SlaterCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
                    double Rkac=ec*ed*cd.Coulomb_R0(la,lc);
                    if (c->kappa>0)
                    {
                        int k2=2*c->kappa+1;
                        Rkac-=(ec+ed)*k2*cd.Coulomb_R0(la,lc-1);
                        Rkac+=k2*k2*cd.Coulomb_R0(la,lc-2);
                    }
                    Jab(ic,id)=Akac(1)*Rkac*norm;
                }
            }
        }
    }
    return J;
};

ERI4 Small_IntegralEngine::MakeDirectSL(const IrrepIEClient* a, const IrrepIEClient* c) const
{
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        double ea=a->es_indices[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            double ec=c->es_indices[ic-1];
            //loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            // Small sector kappa is negative of the large sector kappa, which shifts la but obly for the angular integrals.
            RVec Akac=Coulomb_AngularIntegrals(Omega_kQN::l(-a->kappa),lc,0,0);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                double eb=a->es_indices[ib-1];
                //loop_3(a->es_indices[ib-1]);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << Jab(ic,id) << endl;    
                        assert(false);
                    }
                    double ed=c->es_indices[id-1];
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    SlaterCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
                    double Rkac=ea*eb*cd.Coulomb_R0(la,lc);
                    if (a->kappa>0)
                    {
                        int k2=2*a->kappa+1;
                        Rkac-=(ea+eb)*k2*cd.Coulomb_R0(la-1,lc);
                        Rkac+=k2*k2*cd.Coulomb_R0(la-2,lc);
                    }
                    Jab(ic,id)=Akac(1)*Rkac*norm;
                }
            }
        }
    }
    return J;
}

ERI4 Small_IntegralEngine::MakeDirectSS(const IrrepIEClient* a, const IrrepIEClient* c) const
{
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        double ea=a->es_indices[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            double ec=c->es_indices[ic-1];
            //loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            // Small sector kappa is negative of the large sector kappa, which shifts la &l c but obly for the angular integrals.
            RVec Akac=Coulomb_AngularIntegrals(Omega_kQN::l(-a->kappa),Omega_kQN::l(-c->kappa),0,0);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                double eb=a->es_indices[ib-1];
                //loop_3(a->es_indices[ib-1]);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << Jab(ic,id) << endl;    
                        assert(false);
                    }
                    double ed=c->es_indices[id-1];
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    SlaterCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
                    double Rab0=ea*eb*cd.Coulomb_R0(la,lc);
                    if (a->kappa>0)
                    {
                        int k2=2*a->kappa+1;
                        Rab0-=(ea+eb)*k2*cd.Coulomb_R0(la-1,lc);
                        Rab0+=k2*k2*cd.Coulomb_R0(la-2,lc);
                    }
                    double R0=ec*ed*Rab0;
                    if (c->kappa>0)
                    {
                        double Rab1=ea*eb*cd.Coulomb_R0(la,lc-1);
                        double Rab2=ea*eb*cd.Coulomb_R0(la,lc-2);
                        if (a->kappa>0)
                        {
                            int k2a=2*a->kappa+1;
                            Rab1-=(ea+eb)*k2a*cd.Coulomb_R0(la-1,lc-1);
                            Rab1+=k2a*k2a*cd.Coulomb_R0(la-2,lc-1);
                            Rab2-=(ea+eb)*k2a*cd.Coulomb_R0(la-1,lc-2);
                            Rab2+=k2a*k2a*cd.Coulomb_R0(la-2,lc-2);
                        }
                        int k2c=2*c->kappa+1;
                        R0-=(ec+ed)*k2c*Rab1;
                        R0+=k2c*k2c*Rab2;
                    }
                    Jab(ic,id)=Akac(1)*R0*norm;
                }
            }
        }
    }
    return J;
}


} //namespace
