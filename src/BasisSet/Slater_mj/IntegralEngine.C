// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater_mj/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_mj/IEClient.H" 
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"
#include <iomanip>

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

DiracIntegralEngine::SMat DiracIntegralEngine::merge_off_diag(const SMat& ls)
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
    SMat kl=2.0*itsLargeIE->MakeKinetic(da->itsLargeIEC);
    //SMat ks=itsSmallIE->MakeKinetic(da->itsSmallIEC);
    return merge_off_diag(kl);
}
DiracIntegralEngine::SMat DiracIntegralEngine::MakeNuclear  (iec_t* a, const Cluster& cl) const
{
    auto da=dcast(a);
    SMat vl=itsLargeIE->MakeNuclear(da->itsLargeIEC,cl);
    SMat vs=itsSmallIE->MakeNuclear(da->itsSmallIEC,cl);
    return merge_diag(vl,vs);
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
    assert(false);
    return ERI4();    
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


// Overlap gets called with l*2
double Small_IntegralEngine::Overlap  (double ea , double eb,size_t l2) const
{
    assert(l2%2==0);
    return 2.0*Slater_m::IntegralEngine::Kinetic(ea,eb,l2/2); //Kinetic already has 4*Pi
}

//  This is anew one <a|p^2/r|b> !
double Small_IntegralEngine::Nuclear(double ea, double eb,size_t l) const
{
    double ab=ea+eb;
    int na=l+1,nb=l+1;
    int n=na+nb;
    double Term1= (na*nb+l*(l+1))* SlaterIntegral(ab,n-3); 
    double Term2=-(na*eb+nb*ea)  * SlaterIntegral(ab,n-2);
    double Term3=        ea*eb   * SlaterIntegral(ab,n-1);
    return Term1+Term2+Term3;
}

} //namespace
