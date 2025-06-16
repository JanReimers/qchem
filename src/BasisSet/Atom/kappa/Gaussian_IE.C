// File: Atom/kappa/Gaussian_IE.C  Integral Engine for RKB Gaussians.


#include "kappa/Gaussian_IE.H"
#include "radial/Gaussian/Integrals.H"
// #include "ERI4.H"
// #include <Imp/Symmetry/Okmj.H>
// #include <iostream>

namespace Atom_kappa
{
namespace Gaussian
{


// ERI4 DiracIntegralEngine::merge_diag(const ERI4& LLLL,const ERI4& LLSS,const ERI4& SSLL,const ERI4& SSSS)
// {
//     size_t Nl=LLLL.Nab();
//     size_t Ns=SSSS.Nab();
//     ERI4 J(Nl+Ns,Nl+Ns);
//     for (auto i:LLLL.rows())
//         for (auto j:LLLL.cols(i))
//             J(i,j)=merge_diag(LLLL(i,j),LLSS(i,j));
//     for (auto i:SSSS.rows())
//         for (auto j:SSSS.cols(i))
//             J(Nl+i,Nl+j)=merge_diag(SSLL(i,j),SSSS(i,j));
//     return J;
// }
// ERI4 DiracIntegralEngine::merge_off_diag(const ERI4& LLLL,const M4& LSLS,const M4& SLSL,const ERI4& SSSS)
// {
//     size_t Nl=LLLL.Nab();
//     size_t Ns=SSSS.Nab();
//     ERI4 K(Nl+Ns,Nl+Ns);
//     SMat LLSS(Ns);
//     Fill(LLSS,0.0);
//     for (auto i:LLLL.rows())
//         for (auto j:LLLL.cols(i))
//             K(i,j)=merge_diag(LLLL(i,j),LLSS);
//     for (auto i:SSSS.rows())
//         for (auto j:SSSS.cols(i))
//             K(Nl+i,Ns+j)=merge_diag(LLSS,SSSS(i,j));
//     for (auto i:LSLS.rows())
//         for (auto j:LSLS.cols(i))
//             K(i,Ns+j)=merge_off_diag(LSLS(i,j));   
//     // for (auto i:SLSL.rows())
//     //     for (auto j:SLSL.cols(i))
//     //         K(Nl+i,j)=SLSL(i,j);                  
//     return K;
// }
// ERI4 DiracIntegralEngine::MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const
// {
//     auto da=dcast(a);
//     auto dc=dcast(c);
//     assert(da->itsLargeIEC->kappa==da->itsSmallIEC->kappa);
//     assert(dc->itsLargeIEC->kappa==dc->itsSmallIEC->kappa);
//     ERI4 JLLLL=itsLargeIE->MakeDirect(da->itsLargeIEC,dc->itsLargeIEC);
//     ERI4 JLLSS=itsSmallIE->MakeDirect(da->itsLargeIEC,dc->itsSmallIEC);
//     ERI4 JSSLL=itsSmallIE->MakeDirect(da->itsSmallIEC,dc->itsLargeIEC);
//     ERI4 JSSSS=itsSmallIE->MakeDirect(da->itsSmallIEC,dc->itsSmallIEC);  
//     return merge_diag(JLLLL,JLLSS,JSSLL,JSSSS); 
// }
// ERI4 DiracIntegralEngine::MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const
// {
//     auto da=dcast(a);
//     auto db=dcast(b);
//     auto sie=dynamic_cast<const Small_IntegralEngine*>(itsSmallIE);
//     assert(da->itsLargeIEC->kappa==da->itsSmallIEC->kappa);
//     assert(db->itsLargeIEC->kappa==db->itsSmallIEC->kappa);
//     ERI4 KLLLL=itsLargeIE->MakeExchange(da->itsLargeIEC,db->itsLargeIEC);
//     M4   KLSLS=sie->MakeExchangeLS(da->itsLargeIEC,db->itsSmallIEC);
//     M4   KSLSL=sie->MakeExchangeSL(da->itsSmallIEC,db->itsLargeIEC);
//     ERI4 KSSSS=sie->MakeExchangeSS(da->itsSmallIEC,db->itsSmallIEC);
//     // std::cout << "KLLLL(2,2)=" << KLLLL(2,2) << std::endl;
//     // std::cout << "KLSLS(2,2)=" << KLSLS(2,2) << std::endl;
//     // std::cout << "KSLSL(2,2)=" << KSLSL(2,2) << std::endl;
//     // std::cout << "KSSSS(2,2)=" << KSSSS(2,2) << std::endl;  
//     return merge_off_diag(KLLLL,KLSLS,KSLSL,KSSSS);  
// }
// ERI4 Small_IntegralEngine::MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const
// {
//     auto da=dcast(a);
//     auto dc=dcast(c);  
//     if (da->Large() && !dc->Large())
//         return MakeDirectLS(da,dc);
//     else if (!da->Large() && dc->Large())
//         return MakeDirectSL(da,dc);
//     else if (!da->Large() && !dc->Large())
//         return MakeDirectSS(da,dc);
//     assert(false);
//     return ERI4();
// }
// ERI4 Small_IntegralEngine::MakeDirectLS(const IrrepIEClient* a, const IrrepIEClient* c) const
// {
//     assert(a);
//     assert(c);
//     size_t Na=a->size(), Nc=c->size();
//     ERI4 J(Na,Nc);
//     for (size_t ia:a->indices())
//     {
//         double ea=a->es[ia-1];
//         //loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
//         for (size_t ic:c->indices())
//         {
//             double ec=c->es_indices[ic-1];
//             //loop_2(c->es_indices[ic-1]);
//             int la=a->l, lc=c->l;
//             // Small sector kappa is negative of the large sector kappa, which shifts lc but obly for the angular integrals.
//             RVec Akac=Coulomb_AngularIntegrals(la,Omega_kQN::l(-c->kappa),0,0);
//             for (size_t ib:a->indices())
//             {
//                 if (ib<ia) continue; 
//                 SMat& Jab=J(ia,ib);
//                 double eb=a->es_indices[ib-1];
//                 //loop_3(a->es_indices[ib-1]);
//                 for (size_t id:c->indices())
//                 {
//                     if (id<ic) continue;
//                     if (Jab(ic,id)!=0.0)
//                     {
//                         cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                         cout << Jab(ic,id) << endl;    
//                         assert(false);
//                     }
//                     double ed=c->es[id-1];
//                     double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
//                     Gaussian::RkEngine cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
//                     double Rkac=ec*ed*cd.Coulomb_R0(la,lc);
//                     if (c->kappa>0)
//                     {
//                         int k2=2*c->kappa+1;
//                         Rkac-=(ec+ed)*k2*cd.Coulomb_R0(la,lc-1);
//                         Rkac+=k2*k2*cd.Coulomb_R0(la,lc-2);
//                     }
//                     Jab(ic,id)=Akac(1)*Rkac*norm;
//                 }
//             }
//         }
//     }
//     return J;
// };
// ERI4 Small_IntegralEngine::MakeDirectSL(const IrrepIEClient* a, const IrrepIEClient* c) const
// {
//     assert(a);
//     assert(c);
//     size_t Na=a->size(), Nc=c->size();
//     ERI4 J(Na,Nc);
//     for (size_t ia:a->indices())
//     {
//         double ea=a->es[ia-1];
//         //loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
//         for (size_t ic:c->indices())
//         {
//             double ec=c->es[ic-1];
//             //loop_2(c->es_indices[ic-1]);
//             int la=a->l, lc=c->l;
//             // Small sector kappa is negative of the large sector kappa, which shifts la but obly for the angular integrals.
//             RVec Akac=Coulomb_AngularIntegrals(Omega_kQN::l(-a->kappa),lc,0,0);
//             for (size_t ib:a->indices())
//             {
//                 if (ib<ia) continue; 
//                 SMat& Jab=J(ia,ib);
//                 double eb=a->es[ib-1];
//                 //loop_3(a->es_indices[ib-1]);
//                 for (size_t id:c->indices())
//                 {
//                     if (id<ic) continue;
//                     if (Jab(ic,id)!=0.0)
//                     {
//                         cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                         cout << Jab(ic,id) << endl;    
//                         assert(false);
//                     }
//                     double ed=c->es[id-1];
//                     double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
//                     Gaussian::RkEngine cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
//                     double Rkac=ea*eb*cd.Coulomb_R0(la,lc);
//                     if (a->kappa>0)
//                     {
//                         int k2=2*a->kappa+1;
//                         Rkac-=(ea+eb)*k2*cd.Coulomb_R0(la-1,lc);
//                         Rkac+=k2*k2*cd.Coulomb_R0(la-2,lc);
//                     }
//                     Jab(ic,id)=Akac(1)*Rkac*norm;
//                 }
//             }
//         }
//     }
//     return J;
// }
// ERI4 Small_IntegralEngine::MakeDirectSS(const IrrepIEClient* a, const IrrepIEClient* c) const
// {
//     assert(a);
//     assert(c);
//     size_t Na=a->size(), Nc=c->size();
//     ERI4 J(Na,Nc);
//     for (size_t ia:a->indices())
//     {
//         double ea=a->es[ia-1];
//         //loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
//         for (size_t ic:c->indices())
//         {
//             double ec=c->es[ic-1];
//             //loop_2(c->es_indices[ic-1]);
//             int la=a->l, lc=c->l;
//             // Small sector kappa is negative of the large sector kappa, which shifts la &l c but obly for the angular integrals.
//             RVec Akac=Coulomb_AngularIntegrals(Omega_kQN::l(-a->kappa),Omega_kQN::l(-c->kappa),0,0);
//             for (size_t ib:a->indices())
//             {
//                 if (ib<ia) continue; 
//                 SMat& Jab=J(ia,ib);
//                 double eb=a->es[ib-1];
//                 //loop_3(a->es_indices[ib-1]);
//                 for (size_t id:c->indices())
//                 {
//                     if (id<ic) continue;
//                     if (Jab(ic,id)!=0.0)
//                     {
//                         cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                         cout << Jab(ic,id) << endl;    
//                         assert(false);
//                     }
//                     double ed=c->es[id-1];
//                     double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
//                     Gaussian::RkEngine cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
//                     double Rab0=ea*eb*cd.Coulomb_R0(la,lc);
//                     if (a->kappa>0)
//                     {
//                         int k2=2*a->kappa+1;
//                         Rab0-=(ea+eb)*k2*cd.Coulomb_R0(la-1,lc);
//                         Rab0+=k2*k2*cd.Coulomb_R0(la-2,lc);
//                     }
//                     double R0=ec*ed*Rab0;
//                     if (c->kappa>0)
//                     {
//                         double Rab1=ea*eb*cd.Coulomb_R0(la,lc-1);
//                         double Rab2=ea*eb*cd.Coulomb_R0(la,lc-2);
//                         if (a->kappa>0)
//                         {
//                             int k2a=2*a->kappa+1;
//                             Rab1-=(ea+eb)*k2a*cd.Coulomb_R0(la-1,lc-1);
//                             Rab1+=k2a*k2a*cd.Coulomb_R0(la-2,lc-1);
//                             Rab2-=(ea+eb)*k2a*cd.Coulomb_R0(la-1,lc-2);
//                             Rab2+=k2a*k2a*cd.Coulomb_R0(la-2,lc-2);
//                         }
//                         int k2c=2*c->kappa+1;
//                         R0-=(ec+ed)*k2c*Rab1;
//                         R0+=k2c*k2c*Rab2;
//                     }
//                     Jab(ic,id)=Akac(1)*R0*norm;
//                 }
//             }
//         }
//     }
//     return J;
// }
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
// M4 Small_IntegralEngine::MakeExchangeLS(const IrrepIEClient* a, const IrrepIEClient* c) const
// {
//     assert(a);
//     assert(c);
//     size_t Na=a->size(), Nc=c->size();
//     M4 K(Na,Nc);
//     for (size_t ia:a->indices())
//     {
//         double ea=a->es[ia-1];
//         //loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
//         double na=a->ns(ia);
//         for (size_t ic:c->indices())
//         {
//             int la=a->l, lc=c->l;
//             RVec Akac=ExchangeAngularIntegrals(la,Omega_kQN::l(-c->kappa),0,0);
//             double nac=na*c->ns(ic);
//             for (size_t ib:a->indices())
//             {
//                 Mat& Kab=K(ia,ib);
//                 double eb=a->es[ib-1];
//                 double ec=c->es[ic-1];
//                 // loop_2(a->es_indices[ib-1]);
//                 // loop_3(c->es_indices[ic-1]);
//                 double nacb=nac*a->ns(ib);
//                 for (size_t id:c->indices())
//                 {
//                     double ed=c->es[id-1];
//                     double norm=nacb*c->ns(id);
//                     Gaussian::RkEngine cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
//                     RVec Rkac=ec*ed*cd.ExchangeRk(la,lc);
//                     //RVec RKac=loop_4_exchange(c->es_indices[id-1],la,lc);
//                      if (c->kappa>0)
//                     {
//                         int k2=2*c->kappa+1;
//                         Rkac-=RVec((ec+ed)*k2*cd.ExchangeRk(la,lc-1));
//                         Rkac+=RVec(k2*k2*cd.ExchangeRk(la,lc-2));
//                     }
//                     Kab(ic,id)=Akac*Rkac*norm; //THis whole block is off diagonal.                   
//                 }
//             }
//         }
//     }
//     return K;
// }
// M4 Small_IntegralEngine::MakeExchangeSL(const IrrepIEClient* a, const IrrepIEClient* c) const
// {
//     assert(a);
//     assert(c);
//     size_t Na=a->size(), Nc=c->size();
//     M4 K(Na,Nc);
//     for (size_t ia:a->indices())
//     {
//         double ea=a->es[ia-1];
//         //loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
//         double na=a->ns(ia);
//         for (size_t ic:c->indices())
//         {
//             int la=a->l, lc=c->l;
//             RVec Akac=ExchangeAngularIntegrals(Omega_kQN::l(-a->kappa),lc,0,0);
//             double nac=na*c->ns(ic);
//             for (size_t ib:a->indices())
//             {
//                 Mat& Kab=K(ia,ib);
//                 double eb=a->es[ib-1];
//                 double ec=c->es[ic-1];
//                 // loop_2(a->es_indices[ib-1]);
//                 // loop_3(c->es_indices[ic-1]);
//                 double nacb=nac*a->ns(ib);
//                 for (size_t id:c->indices())
//                 {
//                     double ed=c->es[id-1];
//                     double norm=nacb*c->ns(id);
//                     Gaussian::RkEngine cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
//                     RVec Rkac=ea*eb*cd.ExchangeRk(la,lc);
//                     //RVec RKac=loop_4_exchange(c->es_indices[id-1],la,lc);
//                      if (a->kappa>0)
//                     {
//                         int k2=2*a->kappa+1;
//                         Rkac-=RVec((ea+eb)*k2*cd.ExchangeRk(la-1,lc));
//                         Rkac+=RVec(k2*k2*cd.ExchangeRk(la-2,lc));
//                     }
//                     Kab(ic,id)=Akac*Rkac*norm; //THis whole block is off diagonal.
//                 }
//             }
//         }
//     }
//     return K;
// }
// ERI4 Small_IntegralEngine::MakeExchangeSS(const IrrepIEClient*a, const IrrepIEClient* c) const
// {
//     assert(a);
//     assert(c);
//     size_t Na=a->size(), Nc=c->size();
//     ERI4 K(Na,Nc);
//     for (size_t ia:a->indices())
//     {
//         double ea=a->es[ia-1];
//         //loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
//         double na=a->ns(ia);
//         for (size_t ic:c->indices())
//         {
//             int la=a->l, lc=c->l;
//             int Ala=Omega_kQN::l(-a->kappa), Alc=Omega_kQN::l(-c->kappa);
//             RVec Akac=ExchangeAngularIntegrals(Ala,Alc,0,0);
//             double nac=na*c->ns(ic);
//             for (size_t ib:a->indices(ia))
//             {
//                 SMat& Kab=K(ia,ib);
//                 double eb=a->es[ib-1];
//                 double ec=c->es[ic-1];
//                 // loop_2(a->es_indices[ib-1]);
//                 // loop_3(c->es_indices[ic-1]);
//                 double nacb=nac*a->ns(ib);
//                 for (size_t id:c->indices())
//                 {
//                     double ed=c->es[id-1];
//                     double norm=nacb*c->ns(id);
//                     Gaussian::RkEngine cd(ea+eb,ec+ed,LMax(ia,ib,ic,id));
//                     RVec Rabk=ea*eb*cd.ExchangeRk(la,lc);
//                     if (a->kappa>0)
//                     {
//                         int k2=2*a->kappa+1;
//                         Rabk-=RVec((ea+eb)*k2*cd.ExchangeRk(la-1,lc));
//                         Rabk+=RVec(k2*k2*cd.ExchangeRk(la-2,lc));
//                     }
//                     RVec Rk=ec*ed*Rabk;
//                     if (c->kappa>0)
//                     {
//                         RVec Rab1=ea*eb*cd.ExchangeRk(la,lc-1);
//                         RVec Rab2=ea*eb*cd.ExchangeRk(la,lc-2);
//                         if (a->kappa>0)
//                         {
//                             int k2a=2*a->kappa+1;
//                             Rab1-=RVec((ea+eb)*k2a*cd.ExchangeRk(la-1,lc-1));
//                             Rab1+=RVec(k2a*k2a*cd.ExchangeRk(la-2,lc-1));
//                             Rab2-=RVec((ea+eb)*k2a*cd.ExchangeRk(la-1,lc-2));
//                             Rab2+=RVec(k2a*k2a*cd.ExchangeRk(la-2,lc-2));
//                         }
//                         int k2c=2*c->kappa+1;
//                         Rk-=RVec((ec+ed)*k2c*Rab1);
//                         Rk+=RVec(k2c*k2c*Rab2);
//                     }
//                     // std::cout << "Rk=" << Rk << std::endl;
//                     // std::cout << "Akac=" << Akac << std::endl;
//                     if (ic==id)
//                         Kab(ic,id)=Akac*Rk*norm; 
//                     else if (id<ic)
//                         Kab(id,ic)+=0.5*Akac*Rk*norm; 
//                     else
//                         Kab(ic,id)+=0.5*Akac*Rk*norm;
//                 }
//             }
//         }
//     }
//     return K;
// }
// 

template <class T> double Orbital_RKBS_IE<T>::Inv_r1(double ea , double eb,size_t l_total) const
{
    assert(l_total==0);
    // +2 from the d/dr and -1 from the 1/r potential = l_total+1
    return 4*ea*eb*::Gaussian::Integral(ea+eb,l_total+1); //Don't count the r^2 in dr^3
}

template class Orbital_RKBS_IE<double>;

}} //namespace
