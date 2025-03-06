// File: SphericalGaussian/IntegralEngine.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"
#include <Imp/Symmetry/OkmjQN.H>
#include <iostream>

using std::cout;
using std::endl;

namespace SphericalGaussian_RKB
{

DiracIntegralEngine::DiracIntegralEngine()
    : itsLargeIE(new SphericalGaussian_m::IntegralEngine())
    , itsSmallIE(new Small_IntegralEngine())
{
    assert(itsLargeIE);
    assert(itsSmallIE);
}
void DiracIntegralEngine::Append(const ::IrrepIEClient* iec)
{
    AnalyticIE<double>::Append(iec);
    const Dirac_IrrepIEClient* diec=dynamic_cast<const Dirac_IrrepIEClient*>(iec);
    itsLargeIE->Append(diec->itsLargeIEC);
    itsSmallIE->Append(diec->itsSmallIEC);
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

ERI4 DiracIntegralEngine::merge_diag(const ERI4& LLLL,const ERI4& LLSS,const ERI4& SSLL,const ERI4& SSSS)
{
    size_t Nl=LLLL.Nab();
    size_t Ns=SSSS.Nab();
    ERI4 J(Nl+Ns,Nl+Ns);
    for (auto i:LLLL.rows())
        for (auto j:LLLL.cols(i))
            J(i,j)=merge_diag(LLLL(i,j),LLSS(i,j));
    for (auto i:SSSS.rows())
        for (auto j:SSSS.cols(i))
            J(Nl+i,Nl+j)=merge_diag(SSLL(i,j),SSSS(i,j));
    return J;
}

ERI4 DiracIntegralEngine::merge_off_diag(const ERI4& LLLL,const M4& LSLS,const M4& SLSL,const ERI4& SSSS)
{
    size_t Nl=LLLL.Nab();
    size_t Ns=SSSS.Nab();
    ERI4 K(Nl+Ns,Nl+Ns);
    SMat LLSS(Ns);
    Fill(LLSS,0.0);
    for (auto i:LLLL.rows())
        for (auto j:LLLL.cols(i))
            K(i,j)=merge_diag(LLLL(i,j),LLSS);
    for (auto i:SSSS.rows())
        for (auto j:SSSS.cols(i))
            K(Nl+i,Ns+j)=merge_diag(LLSS,SSSS(i,j));
    for (auto i:LSLS.rows())
        for (auto j:LSLS.cols(i))
            K(i,Ns+j)=merge_off_diag(LSLS(i,j));   
    // for (auto i:SLSL.rows())
    //     for (auto j:SLSL.cols(i))
    //         K(Nl+i,j)=SLSL(i,j);        
            
    return K;
}

DiracIntegralEngine::SMat DiracIntegralEngine::MakeOverlap  (iec_t* a) const
{
    auto da=dcast(a);
    SMat ol=itsLargeIE->MakeOverlap(da->itsLargeIEC);
    SMat os=itsSmallIE->MakeOverlap(da->itsSmallIEC);
    return merge_diag(ol,os);
}

// DiracIntegralEngine::SMat DiracIntegralEngine::MakeKinetic  (iec_t* a) const
// {
//     auto da=dcast(a);
//     Mat kls=-2.0*itsLargeIE->MakeKinetic(da->itsLargeIEC,da->itsSmallIEC);
//     //std::cout << "kls=" << kls << std::endl;
//     return merge_off_diag(kls);
// }

// DiracIntegralEngine::Mat DiracIntegralEngine::MakeKinetic(iec_t* a,iec_t* b) const
// {
//     assert(false);
//     return Mat();
// }

// DiracIntegralEngine::SMat DiracIntegralEngine::MakeNuclear  (iec_t* a, const Cluster& cl) const
// {
//     auto da=dcast(a);
//     SMat vl=itsLargeIE->MakeNuclear(da->itsLargeIEC,cl);
//     SMat vs=itsSmallIE->MakeNuclear(da->itsSmallIEC,cl);
//     return merge_diag(vl,vs);
// }
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
    assert(da->itsLargeIEC->kappa==da->itsSmallIEC->kappa);
    assert(dc->itsLargeIEC->kappa==dc->itsSmallIEC->kappa);
    ERI4 JLLLL=itsLargeIE->MakeDirect(da->itsLargeIEC,dc->itsLargeIEC);
    ERI4 JLLSS=itsSmallIE->MakeDirect(da->itsLargeIEC,dc->itsSmallIEC);
    ERI4 JSSLL=itsSmallIE->MakeDirect(da->itsSmallIEC,dc->itsLargeIEC);
    ERI4 JSSSS=itsSmallIE->MakeDirect(da->itsSmallIEC,dc->itsSmallIEC);
    
    return merge_diag(JLLLL,JLLSS,JSSLL,JSSSS); 
}
ERI4 DiracIntegralEngine::MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const
{
    auto da=dcast(a);
    auto db=dcast(b);
    auto sie=dynamic_cast<const Small_IntegralEngine*>(itsSmallIE);
    assert(da->itsLargeIEC->kappa==da->itsSmallIEC->kappa);
    assert(db->itsLargeIEC->kappa==db->itsSmallIEC->kappa);
    ERI4 KLLLL=itsLargeIE->MakeExchange(da->itsLargeIEC,db->itsLargeIEC);
    M4   KLSLS=sie->MakeExchangeLS(da->itsLargeIEC,db->itsSmallIEC);
    M4   KSLSL=sie->MakeExchangeSL(da->itsSmallIEC,db->itsLargeIEC);
    ERI4 KSSSS=sie->MakeExchangeSS(da->itsSmallIEC,db->itsSmallIEC);

    // std::cout << "KLLLL(2,2)=" << KLLLL(2,2) << std::endl;
    // std::cout << "KLSLS(2,2)=" << KLSLS(2,2) << std::endl;
    // std::cout << "KSLSL(2,2)=" << KSLSL(2,2) << std::endl;
    // std::cout << "KSSSS(2,2)=" << KSSSS(2,2) << std::endl;
    
    return merge_off_diag(KLLLL,KLSLS,KSLSL,KSSSS);  
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
    return 2.0*SphericalGaussian::IntegralEngine::Kinetic(ea,eb,l2/2); //Kinetic already has 4*Pi
}

//  This is anew one <a|p^2/r|b> !
//  For this one we actually need to know the sign kappa = l or -l-1
double Small_IntegralEngine::Nuclear(double ea, double eb,size_t l) const
{
    assert(l==0);
    //int kappa = -l -1;
    return 4*ea*eb*GaussianIntegral(ea+eb,l+1); //Don't count the r^2 in dr^3
   
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
        double ea=a->es[ia-1];
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
                    double ed=c->es[id-1];
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    SphericalGaussianCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
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
        double ea=a->es[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            double ec=c->es[ic-1];
            //loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            // Small sector kappa is negative of the large sector kappa, which shifts la but obly for the angular integrals.
            RVec Akac=Coulomb_AngularIntegrals(Omega_kQN::l(-a->kappa),lc,0,0);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                double eb=a->es[ib-1];
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
                    double ed=c->es[id-1];
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    SphericalGaussianCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
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
        double ea=a->es[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            double ec=c->es[ic-1];
            //loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            // Small sector kappa is negative of the large sector kappa, which shifts la &l c but obly for the angular integrals.
            RVec Akac=Coulomb_AngularIntegrals(Omega_kQN::l(-a->kappa),Omega_kQN::l(-c->kappa),0,0);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                double eb=a->es[ib-1];
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
                    double ed=c->es[id-1];
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    SphericalGaussianCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
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

// ERI4 Small_IntegralEngine::MakeExchange  (const ::IrrepIEClient* a, const ::IrrepIEClient* b) const
// {
//     auto da=dcast(a);
//     auto db=dcast(b);  
//     if (da->Large() && !db->Large())
//         return MakeExchangeLS(da,db);
//     else if (!da->Large() && db->Large())
//         return MakeExchangeSL(da,db);
//     else if (!da->Large() && !db->Large())
//         return MakeExchangeSS(da,db);
//     assert(false);
//     return ERI4();
// }

M4 Small_IntegralEngine::MakeExchangeLS(const IrrepIEClient* a, const IrrepIEClient* c) const
{
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    M4 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        double ea=a->es[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        double na=a->ns(ia);
        for (size_t ic:c->indices())
        {
            int la=a->l, lc=c->l;
            RVec Akac=ExchangeAngularIntegrals(la,Omega_kQN::l(-c->kappa),0,0);
            double nac=na*c->ns(ic);
            for (size_t ib:a->indices())
            {
                Mat& Kab=K(ia,ib);
                double eb=a->es[ib-1];
                double ec=c->es[ic-1];
                // loop_2(a->es_indices[ib-1]);
                // loop_3(c->es_indices[ic-1]);
                double nacb=nac*a->ns(ib);
                for (size_t id:c->indices())
                {
                    double ed=c->es[id-1];
                    double norm=nacb*c->ns(id);
                    SphericalGaussianCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
                    RVec Rkac=ec*ed*cd.ExchangeRk(la,lc);
                    //RVec RKac=loop_4_exchange(c->es_indices[id-1],la,lc);
                     if (c->kappa>0)
                    {
                        int k2=2*c->kappa+1;
                        Rkac-=RVec((ec+ed)*k2*cd.ExchangeRk(la,lc-1));
                        Rkac+=RVec(k2*k2*cd.ExchangeRk(la,lc-2));
                    }
                    Kab(ic,id)=Akac*Rkac*norm; //THis whole block is off diagonal.
                    
                }
            }
        }
    }

    return K;
}

M4 Small_IntegralEngine::MakeExchangeSL(const IrrepIEClient* a, const IrrepIEClient* c) const
{
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    M4 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        double ea=a->es[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        double na=a->ns(ia);
        for (size_t ic:c->indices())
        {
            int la=a->l, lc=c->l;
            RVec Akac=ExchangeAngularIntegrals(Omega_kQN::l(-a->kappa),lc,0,0);
            double nac=na*c->ns(ic);
            for (size_t ib:a->indices())
            {
                Mat& Kab=K(ia,ib);
                double eb=a->es[ib-1];
                double ec=c->es[ic-1];
                // loop_2(a->es_indices[ib-1]);
                // loop_3(c->es_indices[ic-1]);
                double nacb=nac*a->ns(ib);
                for (size_t id:c->indices())
                {
                    double ed=c->es[id-1];
                    double norm=nacb*c->ns(id);
                    SphericalGaussianCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
                    RVec Rkac=ea*eb*cd.ExchangeRk(la,lc);
                    //RVec RKac=loop_4_exchange(c->es_indices[id-1],la,lc);
                     if (a->kappa>0)
                    {
                        int k2=2*a->kappa+1;
                        Rkac-=RVec((ea+eb)*k2*cd.ExchangeRk(la-1,lc));
                        Rkac+=RVec(k2*k2*cd.ExchangeRk(la-2,lc));
                    }
                    Kab(ic,id)=Akac*Rkac*norm; //THis whole block is off diagonal.
                }
            }
        }
    }
    return K;
}

ERI4 Small_IntegralEngine::MakeExchangeSS(const IrrepIEClient*a, const IrrepIEClient* c) const
{
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        double ea=a->es[ia-1];
        //loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        double na=a->ns(ia);
        for (size_t ic:c->indices())
        {
            int la=a->l, lc=c->l;
            int Ala=Omega_kQN::l(-a->kappa), Alc=Omega_kQN::l(-c->kappa);
            RVec Akac=ExchangeAngularIntegrals(Ala,Alc,0,0);
            double nac=na*c->ns(ic);
            for (size_t ib:a->indices(ia))
            {
                SMat& Kab=K(ia,ib);
                double eb=a->es[ib-1];
                double ec=c->es[ic-1];
                // loop_2(a->es_indices[ib-1]);
                // loop_3(c->es_indices[ic-1]);
                double nacb=nac*a->ns(ib);
                for (size_t id:c->indices())
                {
                    double ed=c->es[id-1];
                    double norm=nacb*c->ns(id);
                    SphericalGaussianCD cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
                    RVec Rabk=ea*eb*cd.ExchangeRk(la,lc);
                    if (a->kappa>0)
                    {
                        int k2=2*a->kappa+1;
                        Rabk-=RVec((ea+eb)*k2*cd.ExchangeRk(la-1,lc));
                        Rabk+=RVec(k2*k2*cd.ExchangeRk(la-2,lc));
                    }
                    RVec Rk=ec*ed*Rabk;
                    if (c->kappa>0)
                    {
                        RVec Rab1=ea*eb*cd.ExchangeRk(la,lc-1);
                        RVec Rab2=ea*eb*cd.ExchangeRk(la,lc-2);
                        if (a->kappa>0)
                        {
                            int k2a=2*a->kappa+1;
                            Rab1-=RVec((ea+eb)*k2a*cd.ExchangeRk(la-1,lc-1));
                            Rab1+=RVec(k2a*k2a*cd.ExchangeRk(la-2,lc-1));
                            Rab2-=RVec((ea+eb)*k2a*cd.ExchangeRk(la-1,lc-2));
                            Rab2+=RVec(k2a*k2a*cd.ExchangeRk(la-2,lc-2));
                        }
                        int k2c=2*c->kappa+1;
                        Rk-=RVec((ec+ed)*k2c*Rab1);
                        Rk+=RVec(k2c*k2c*Rab2);
                    }
                    // std::cout << "Rk=" << Rk << std::endl;
                    // std::cout << "Akac=" << Akac << std::endl;
                    if (ic==id)
                        Kab(ic,id)=Akac*Rk*norm; 
                    else if (id<ic)
                        Kab(id,ic)+=0.5*Akac*Rk*norm; 
                    else
                        Kab(ic,id)+=0.5*Akac*Rk*norm;
                }
            }
        }
    }
    return K;
}

double Small_IntegralEngine1::Integral(qchem::IType t,double ea , double eb,size_t l) const
{
    if (t==qchem::Overlap1 || t==qchem::Kinetic1)
    {
        assert(l==0);
        double eab=ea+eb;
        size_t l1=l+1;
        return 1.0*(
                (l1*l1 + l*l1) * GaussianIntegral(eab,2*l-2)
                -2*l1 * eab      * GaussianIntegral(eab,2*l  )
                +4*ea*eb       * GaussianIntegral(eab,2*l+2)
            );
    }
    else if (t==qchem::Nuclear1)
    {
        assert(l==0);
        //int kappa = -l -1;
        return 4*ea*eb*GaussianIntegral(ea+eb,l+1); //Don't count the r^2 in dr^3
    }
    assert(false);
    return 0.0;
}


} //namespace
