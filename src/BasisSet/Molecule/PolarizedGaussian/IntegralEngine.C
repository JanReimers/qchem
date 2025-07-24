// File: IntegralEngine.C  Here is where all the integral get calculated.


#include <cassert>
#include <memory>
#include <cmath>
#include "PolarizedGaussian/IEClient.H"
#include "PolarizedGaussian/IntegralEngine.H"
import qchem.Fit_IBS;
import qchem.Irrep_BS;
import qchem.BasisSet.Internal.ERI4;

namespace PolarizedGaussian
{

Fit_IE::Vec Fit_IE::MakeCharge() const
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient*>(this);
    assert(a);
    Vec c(a->size());
    int i=0;
    for (auto r:a->radials)
    {
        c(i+1)=r->GetCharge(a->pols[i])*a->ns(i+1); 
        i++;       
    }

    // 1 line copilot version
    //for (auto i:ab->ns.indices())  c(i)=ab->radials[i-1]->GetCharge(ab->pols[i-1])*ab->ns(i);
    return c;
}

Fit_IE::Mat Fit_IE::MakeRepulsion(const Fit_IBS& _b) const
{   
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient*>(this);
    const IrrepIEClient* b=dynamic_cast<const IrrepIEClient*>(&_b);
    assert(a);
    assert(b); 
    int Na=a->size(),Nb=b->size();
    Mat s(Na,Nb);
    for (size_t ia=0;ia<Na;ia++)
        for (size_t ib=0;ib<Nb;ib++)
            s(ia+1,ib+1)=a->radials[ia]->Integrate(qchem::Repulsion2C,
                b->radials[ib],a->pols[ia],b->pols[ib],cache)*a->ns(ia+1)*b->ns(ib+1);
    assert(!isnan(s));
    return s;
}

SMatrix<double> IE_Common::MakeIntegrals(qchem::IType2C t2C,const Cluster* cl) const
{
    const IrrepIEClient* ab=dynamic_cast<const IrrepIEClient*>(this);
    assert(ab);
    int N=ab->size();
    SMatrix<double> s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=ab->radials[ia]->Integrate(t2C,ab->radials[ib],ab->pols[ia],ab->pols[ib],cache,cl)*ab->ns(ia+1)*ab->ns(ib+1);

    return s;
}

ERI3<double> Orbital_IE::MakeOverlap3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const IrrepIEClient*>(&_c);
    int Nc=c->size();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        SMatrix<double> s=Integrate(qchem::Overlap3C,c->radials[ic],c->pols[ic]);
        s*=c->ns(ic+1);
        s3.push_back(s);
    } 
    return s3;   
}
ERI3<double> Orbital_IE::MakeRepulsion3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const IrrepIEClient*>(&_c);
    int Nc=c->size();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        SMatrix<double> s=Integrate(qchem::Repulsion3C,c->radials[ic],c->pols[ic]);
        s*=c->ns(ic+1);
        s3.push_back(s);
    }    
    return s3;
}
SMatrix<double> Orbital_IE::Integrate(qchem::IType3C type , const RadialFunction* rc, const Polarization& pc) const
{
    auto ab=dynamic_cast<const IrrepIEClient*>(this);
    int N=ab->size();
    SMatrix<double> s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(type,ab->radials[ia],ab->radials[ib],ab->pols[ia],ab->pols[ib],pc,cache)*ab->ns(ia+1)*ab->ns(ib+1);
        
    return s;    
}



ERI4 Orbital_IE::MakeDirect  (const obs_t& _c) const
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(this);
    const IrrepIEClient* c=dynamic_cast<const IrrepIEClient* >(&_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    
    for (size_t ia:a->ns.indices())
        for (size_t ib:a->ns.indices(ia))
        {
            SMatrix<double>& Jab=J(ia,ib);
            for (size_t ic:c->ns.indices())
                for (size_t id:c->ns.indices(ic))
                {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                        assert(c->radials[id-1]);
                        Jab(ic,id)=norm * c->radials[id-1]->Integrate(a->radials[ia-1],a->radials[ib-1],c->radials[ic-1],a->pols[ia-1],a->pols[ib-1],c->pols[ic-1],c->pols[id-1],cache);
                }
        }
    return J;
}

ERI4 Orbital_IE::MakeExchange(const obs_t& _b) const
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(this);
    const IrrepIEClient* b=dynamic_cast<const IrrepIEClient* >(&_b);
    assert(a);
    assert(b);
    size_t Na=a->size(), Nb=b->size();
    ERI4 K(Na,Nb);
    for (size_t ia:a->ns.indices())
        for (size_t ib:b->ns.indices())
           
            for (size_t ic:a->ns.indices(ia))
            {
                SMatrix<double>& Kac=K(ia,ic);
                for (size_t id:b->ns.indices())
                {
                  //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                    double norm=a->ns(ia)*b->ns(ib)*a->ns(ic)*b->ns(id);
                    assert(b->radials[id-1]);
                    if (ib==id)
                        Kac(ib,id)=norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                }        
            }
    return K;
}

} //namespace PolarizedGaussian
