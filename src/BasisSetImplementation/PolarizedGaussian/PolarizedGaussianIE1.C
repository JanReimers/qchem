// File: PolarizedGaussianIE1.C  Here is where all the integral get calculated.


#include "BasisSetImplementation/PolarizedGaussian/IEClient.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianIE1.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Misc/ERI4.H"

//-----------------------------------------------------------------
//
//  Construction zone.
//
PolarizedGaussianIE1::PolarizedGaussianIE1(const blocks_t& _blocks)
    : blocks(_blocks)
    , ns()
    , ons()
{
    for (auto bl:blocks)
        for (auto p:bl->itsPols)
        {
            radials.push_back(bl->itsRadial);
            pols.push_back(p);
        }
   
    size_t N=size();
    ns.SetLimits(N);
    for (size_t i=0;i<N;i++)
        ns(i+1)=radials[i]->Integrate(RadialFunction::Overlap2C,radials[i],pols[i],pols[i],cache);
    ns=1.0/sqrt(ns);
    ons=OuterProduct(ns);
};

PolarizedGaussianIE1::RVec PolarizedGaussianIE1::MakeNormalization() const
{
    return ns;
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
    const PolarizedGaussianIE1* b=dynamic_cast<const PolarizedGaussianIE1*>(ieb);
    assert(b);
    int Na=size(),Nb=b->size();
    Mat s(Na,Nb);
    for (index_t ia=0;ia<Na;ia++)
        for (index_t ib=0;ib<Nb;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Repulsion2C,
                b->radials[ib],pols[ia],b->pols[ib],cache)*ns(ia+1)*b->ns(ib+1);
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

PolarizedGaussianIE1::jk_t PolarizedGaussianIE1::Make4C(const iecv_t& ies) const
{
#ifdef DEBUG
    assert(ies.size()==1);  //Don't have multiple Irrep basis sets for Molecules yet.
    iec_t* ie=ies[0];
    const PolarizedGaussianIEClient* other=dynamic_cast<const PolarizedGaussianIEClient*>(ie);
//    assert(other);
//    assert(other==this);   //Again don't have multiple Irrep basis sets for Molecules yet.
#endif
    
    int N=size();
    ERI4 J(N,-1.0),K;
    std::cout << N << " " << J.itsData.size() <<" " << K.itsData.size() << std::endl;

    for (index_t ia:ns.indices())
        for (index_t ib:ns.indices(ia))
            for (index_t ic:ns.indices())
                for (index_t id:ns.indices(ic))
                {
                    if (J(ia,ib,ic,id)==-1.0)
                    {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=ns(ia)*ns(ib)*ns(ic)*ns(id);
                        assert(radials[id-1]);
                        J(ia,ib,ic,id)=norm * radials[id-1]->Integrate(radials[ia-1],radials[ib-1],radials[ic-1],pols[ia-1],pols[ib-1],pols[ic-1],pols[id-1],cache);
                    }
                }
    return std::make_pair(J,K);
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

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::Integrate(RadialFunction::Types3C type , const RadialFunction* rc, const Polarization& pc) const
{
    int N=size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(type,radials[ia],radials[ib],pols[ia],pols[ib],pc,cache);
        
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

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::Integrate(RadialFunction::Types2C type , const Cluster* cl) const
{
    int N=size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(type,radials[ib],pols[ia],pols[ib],cache,cl);

    Normalize(s);
    return s;
}



std::ostream& PolarizedGaussianIE1::Write(std::ostream& os) const
{
    return os;
}
std::istream& PolarizedGaussianIE1::Read (std::istream& is)
{
    return is;
}



