// File: PolarizedGaussianIE1.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/PolarizedGaussian/IEClient.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Imp/Containers/ERI4.H"

//-----------------------------------------------------------------
//
//  Construction zone.
//
PolarizedGaussianIE1::PolarizedGaussianIE1(const PolarizedGaussianIEClient* iec)
: ns(iec->ns)
, ons(OuterProduct(ns))
{}

PolarizedGaussianIE1::RVec PolarizedGaussianIE1::MakeNormalization(iec_t* iea) const
{
    const PolarizedGaussianIEClient* a=dynamic_cast<const PolarizedGaussianIEClient*>(iea);
    assert(a);
    RVec n(a->size());
    int i=0;
    for (auto r:a->radials)
    {
        n(i+1)=r->Integrate(RadialFunction::Overlap2C,r,a->pols[i],a->pols[i],cache);
        i++;       
    }
    n=1.0/sqrt(n);
    assert(!isnan(n));
    return n;
    
}

PolarizedGaussianIE1::RVec PolarizedGaussianIE1::MakeCharge(iec_t* iea) const
{
    const PolarizedGaussianIEClient* a=dynamic_cast<const PolarizedGaussianIEClient*>(iea);
    assert(a);
    RVec c(a->size());
    int i=0;
    for (auto r:a->radials)
    {
        c(i+1)=r->GetCharge(a->pols[i])*a->ns(i+1); 
        i++;       
    }
    assert(!isnan(c));
    return c;
}


void PolarizedGaussianIE1::Normalize(SMat& s) const
{
    s=DirectMultiply(s,ons);
}

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
AnalyticIE<double>* PolarizedGaussianIE1::Clone() const
{
    return new PolarizedGaussianIE1(*this);
}

//----------------------------------------------------------------------------------------
//
//  2 Center type integrals
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeOverlap(iec_t* a) const
{
    return Integrate(RadialFunction::Overlap2C,a);
}

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeRepulsion(iec_t* iea ) const
{
    return Integrate(RadialFunction::Repulsion2C,iea);
}

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeKinetic(iec_t* a) const
{
    return Integrate(RadialFunction::Kinetic,a);
}
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeNuclear(iec_t* a,const Cluster& cl) const
{
    return Integrate(RadialFunction::Nuclear,a,&cl);
}

PolarizedGaussianIE1::Mat PolarizedGaussianIE1::MakeRepulsion(const IE* iea,const IE* ieb) const
{    
    const PolarizedGaussianIEClient* a=dynamic_cast<const PolarizedGaussianIEClient*>(iea);
    assert(a);
    const PolarizedGaussianIEClient* b=dynamic_cast<const PolarizedGaussianIEClient*>(ieb);
    assert(b);
    int Na=a->size(),Nb=b->size();
    Mat s(Na,Nb);
    for (index_t ia=0;ia<Na;ia++)
        for (index_t ib=0;ib<Nb;ib++)
            s(ia+1,ib+1)=a->radials[ia]->Integrate(RadialFunction::Repulsion2C,
                b->radials[ib],a->pols[ia],b->pols[ib],cache)*a->ns(ia+1)*b->ns(ib+1);
    assert(!isnan(s));
    return s;
}

PolarizedGaussianIE1::Mat PolarizedGaussianIE1::MakeRepulsion(iec_t* iea,iec_t* ieb) const
{    
    const PolarizedGaussianIEClient* a=dynamic_cast<const PolarizedGaussianIEClient*>(iea);;
    assert(a);
    const PolarizedGaussianIEClient* b=dynamic_cast<const PolarizedGaussianIEClient*>(ieb);;
    assert(b);
    int Na=a->size(),Nb=b->size();
    Mat s(Na,Nb);
    for (index_t ia=0;ia<Na;ia++)
        for (index_t ib=0;ib<Nb;ib++)
            s(ia+1,ib+1)=a->radials[ia]->Integrate(RadialFunction::Repulsion2C,
                b->radials[ib],a->pols[ia],b->pols[ib],cache)*a->ns(ia+1)*b->ns(ib+1);
    assert(!isnan(s));
    return s;
}



//------------------------------------------------------------------------------
//
//  3 Centre integrals.
//

PolarizedGaussianIE1::ERI3 PolarizedGaussianIE1::MakeOverlap3C(iec_t* ieab,iec_t* iec) const
{
    const PolarizedGaussianIEClient* c=dynamic_cast<const PolarizedGaussianIEClient*>(iec);
    assert(c);
    int Nc=c->size();
    ERI3 s3;
    for (index_t ic=0;ic<Nc;ic++)
    {
        SMat s=Integrate(RadialFunction::Overlap3C,ieab,c->radials[ic],c->pols[ic]);
        s*=c->ns(ic+1);
        s3.push_back(s);
    } 
    return s3;   
}


PolarizedGaussianIE1::ERI3 PolarizedGaussianIE1::MakeRepulsion3C(iec_t* ieab,iec_t* iec) const
{
    const PolarizedGaussianIEClient* c=dynamic_cast<const PolarizedGaussianIEClient*>(iec);;
    int Nc=c->size();
    ERI3 s3;
    for (index_t ic=0;ic<Nc;ic++)
    {
        SMat s=Integrate(RadialFunction::Repulsion3C,ieab,c->radials[ic],c->pols[ic]);
        s*=c->ns(ic+1);
        s3.push_back(s);
    }    
    return s3;
}

//-----------------------------------------------------------------------------------
//
//  4 centre integrals.
//
PolarizedGaussianIE1::PGparams::PGparams(const iecv_t& iecs)
{
    size_t N=0;
    for (auto iec:iecs) N+=iec->size();
    ns.SetLimits(N);
    for (auto iec:iecs)
    {
        int j=1;
        const PolarizedGaussianIEClient* pg=dynamic_cast<const PolarizedGaussianIEClient*>(iec);
        assert(pg);
        for (size_t i=0;i<pg->size();i++)
        {
            radials.push_back(pg->radials[i]);
            pols.push_back(pg->pols[i]);
            ns(j)=pg->ns(i+1);
            j++;
        }
    }
}



void PolarizedGaussianIE1::Make4C(ERI4& J, ERI4& K,const iecv_t& iecv) const
{
    PGparams abcd(iecv);
    
    int N=abcd.size();
    J.SetSize(N,-1);
    std::cout << N << " " << J.itsData.size() <<" " << K.itsData.size() << std::endl;

    for (index_t ia:abcd.ns.indices())
        for (index_t ib:abcd.ns.indices(ia))
            for (index_t ic:abcd.ns.indices())
                for (index_t id:abcd.ns.indices(ic))
                {
                    if (J(ia,ib,ic,id)==-1.0)
                    {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=abcd.ns(ia)*abcd.ns(ib)*abcd.ns(ic)*abcd.ns(id);
                        assert(abcd.radials[id-1]);
                        J(ia,ib,ic,id)=norm * abcd.radials[id-1]->Integrate(abcd.radials[ia-1],abcd.radials[ib-1],abcd.radials[ic-1],abcd.pols[ia-1],abcd.pols[ib-1],abcd.pols[ic-1],abcd.pols[id-1],cache);
                    }
                }
}


//-------------------------------------------------------------------------
//
//  Internal function for doing most of the double loops.
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::Integrate(RadialFunction::Types3C type ,iec_t* ieab, const RadialFunction* rc, const Polarization& pc) const
{
    const PolarizedGaussianIEClient* ab=dynamic_cast<const PolarizedGaussianIEClient*>(ieab);;
    assert(ab);
    int N=ab->size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(type,ab->radials[ia],ab->radials[ib],ab->pols[ia],ab->pols[ib],pc,cache);
        
    Normalize(s);
    return s;    
}


PolarizedGaussianIE1::SMat PolarizedGaussianIE1::Integrate(RadialFunction::Types2C type ,iec_t* ieab,  const Cluster* cl) const
{
    const PolarizedGaussianIEClient* ab=dynamic_cast<const PolarizedGaussianIEClient*>(ieab);;
    assert(ab);
    int N=ab->size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=ab->radials[ia]->Integrate(type,ab->radials[ib],ab->pols[ia],ab->pols[ib],cache,cl);

    Normalize(s);
    return s;
}



