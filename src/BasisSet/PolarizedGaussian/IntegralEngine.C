// File: IntegralEngine.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/PolarizedGaussian/IEClient.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Imp/Containers/ERI4.H"

namespace PolarizedGaussian
{

const IrrepIEClient* IntegralEngine::dcast(iec_t* iea)
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient*>(iea);
    assert(a);
    return a;
}

IntegralEngine::RVec IntegralEngine::MakeNormalization(iec_t* iea) const
{
    auto a=dcast(iea);
    RVec n(a->size());
    int i=0;
    for (auto r:a->radials)
    {
        n(i+1)=r->Integrate(qchem::Overlap2C,r,a->pols[i],a->pols[i],cache);
        i++;       
    }
    n=1.0/sqrt(n);
    assert(!isnan(n));
    return n;
    
}

IntegralEngine::RVec IntegralEngine::MakeCharge(iec_t* iea) const
{
    auto a=dcast(iea);
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



//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
AnalyticIE<double>* IntegralEngine::Clone() const
{
    return new IntegralEngine(*this);
}

//----------------------------------------------------------------------------------------
//
//  2 Center type integrals
//
IntegralEngine::SMat IntegralEngine::MakeOverlap(iec_t* a) const
{
    return Integrate(qchem::Overlap2C,a);
}

IntegralEngine::SMat IntegralEngine::MakeRepulsion(iec_t* iea ) const
{
    return Integrate(qchem::Repulsion2C,iea);
}

IntegralEngine::SMat IntegralEngine::MakeKinetic(iec_t* a) const
{
    return Integrate(qchem::Kinetic,a);
}
//
IntegralEngine::SMat IntegralEngine::MakeNuclear(iec_t* a,const Cluster& cl) const
{
    return Integrate(qchem::Nuclear,a,&cl);
}

IntegralEngine::Mat IntegralEngine::MakeRepulsion(iec_t* iea,iec_t* ieb) const
{    
    auto a=dcast(iea);
    auto b=dcast(ieb);
    int Na=a->size(),Nb=b->size();
    Mat s(Na,Nb);
    for (index_t ia=0;ia<Na;ia++)
        for (index_t ib=0;ib<Nb;ib++)
            s(ia+1,ib+1)=a->radials[ia]->Integrate(qchem::Repulsion2C,
                b->radials[ib],a->pols[ia],b->pols[ib],cache)*a->ns(ia+1)*b->ns(ib+1);
    assert(!isnan(s));
    return s;
}



//------------------------------------------------------------------------------
//
//  3 Centre integrals.
//

IntegralEngine::ERI3 IntegralEngine::MakeOverlap3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);
    int Nc=c->size();
    ERI3 s3;
    for (index_t ic=0;ic<Nc;ic++)
    {
        SMat s=Integrate(qchem::Overlap3C,ieab,c->radials[ic],c->pols[ic]);
        s*=c->ns(ic+1);
        s3.push_back(s);
    } 
    return s3;   
}


IntegralEngine::ERI3 IntegralEngine::MakeRepulsion3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);
    int Nc=c->size();
    ERI3 s3;
    for (index_t ic=0;ic<Nc;ic++)
    {
        SMat s=Integrate(qchem::Repulsion3C,ieab,c->radials[ic],c->pols[ic]);
        s*=c->ns(ic+1);
        s3.push_back(s);
    }    
    return s3;
}


//-----------------------------------------------------------------------------------
//
//  4 centre integrals.
//

ERIJ IntegralEngine::MakeDirect  (const ::IrrepIEClient* _a, const ::IrrepIEClient* _c,const ::IEClient*) const
{
   // const IEClient* pgiec= dynamic_cast<const IEClient*>(iec);
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(_a);
    const IrrepIEClient* c=dynamic_cast<const IrrepIEClient* >(_c);
    //assert(pgiec);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERIJ J(Na,Nc);
    
    for (index_t ia:a->ns.indices())
        for (index_t ib:a->ns.indices(ia))
            for (index_t ic:c->ns.indices())
                for (index_t id:c->ns.indices(ic))
                {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                        assert(c->radials[id-1]);
                        J(ia,ib,ic,id)=norm * c->radials[id-1]->Integrate(a->radials[ia-1],a->radials[ib-1],c->radials[ic-1],a->pols[ia-1],a->pols[ib-1],c->pols[ic-1],c->pols[id-1],cache);
                }
    return J;
}

ERIK IntegralEngine::MakeExchange(const ::IrrepIEClient* _a, const ::IrrepIEClient* _b,const ::IEClient*) const
{
//    const IEClient* pgiec= dynamic_cast<const IEClient*>(iec);
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(_a);
    const IrrepIEClient* b=dynamic_cast<const IrrepIEClient* >(_b);
//    assert(pgiec);
    assert(a);
    assert(b);
    size_t Na=a->size(), Nb=b->size();
    ERIK K(Na,Nb);
    for (index_t ia:a->ns.indices())
        for (index_t ib:b->ns.indices())
            for (index_t ic:a->ns.indices())
                for (index_t id:b->ns.indices())
                {
                  //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                    double norm=a->ns(ia)*b->ns(ib)*a->ns(ic)*b->ns(id);
                    assert(b->radials[id-1]);
                    K(ia,ic,ib,id)=norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                }        
    return K;
}

//-------------------------------------------------------------------------
//
//  Internal function for doing most of the double loops.
//
IntegralEngine::SMat IntegralEngine::Integrate(qchem::IType3C type ,iec_t* ieab, const RadialFunction* rc, const Polarization& pc) const
{
    auto ab=dcast(ieab);
    int N=ab->size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(type,ab->radials[ia],ab->radials[ib],ab->pols[ia],ab->pols[ib],pc,cache)*ab->ns(ia+1)*ab->ns(ib+1);
        
    return s;    
}


IntegralEngine::SMat IntegralEngine::Integrate(qchem::IType2C type ,iec_t* ieab,  const Cluster* cl) const
{
    auto ab=dcast(ieab);
    int N=ab->size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=ab->radials[ia]->Integrate(type,ab->radials[ib],ab->pols[ia],ab->pols[ib],cache,cl)*ab->ns(ia+1)*ab->ns(ib+1);

    return s;
}

void IntegralEngine::Report(std::ostream& os) const
{
    os << "Polarized Gaussian integral engine cache:" << std::endl;
    cache.Report(os);
}

} //namespace PolarizedGaussian
