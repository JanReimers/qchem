// File: BasisSet/Molecule/PolarizedGaussian/Internal/Imp/IntegralEngine.C  Here is where all the integral2 get calculated.
module;

#include <cassert>
#include <memory>
#include <cmath>
#include "blaze/Math.h"
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IntegralEngine;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Fit_IBS;
import qchem.IrrepBasisSet;
import qchem.BasisSet.Internal.ERI4;


namespace PolarizedGaussian
{

rvec_t Fit_IE::MakeCharge() const
{
    const PGData* a=dynamic_cast<const PGData*>(this);
    assert(a);
    rvec_t c(a->size1());
    int i=0;
    for (auto r:a->radials1)
    {
        c[i]=r->GetCharge(a->pols1[i])*a->ns1[i]; 
        i++;       
    }

    // 1 line copilot version
    //for (auto i:ab->ns.indices())  c(i)=ab->radials[i-1]->GetCharge(ab->pols[i-1])*ab->ns(i);
    return c;
}

rmat_t Fit_IE::MakeRepulsion(const Fit_IBS& _b) const
{   
    const PGData* a=dynamic_cast<const PGData*>(this);
    const PGData* b=dynamic_cast<const PGData*>(&_b);
    assert(a);
    assert(b); 
    int Na=a->size1(),Nb=b->size1();
    rmat_t s(Na,Nb);
    for (size_t ia=0;ia<Na;ia++)
        for (size_t ib=0;ib<Nb;ib++)
            s(ia,ib)=a->radials1[ia]->Integrate(Repulsion2C,
                b->radials1[ib],a->pols1[ia],b->pols1[ib],cache)*a->ns1[ia]*b->ns1[ib];
    assert(!isnan(s));
    return s;
}

rsmat_t IE_Common::MakeIntegrals(IType t2C,const Cluster* cl) const
{
    const PGData* ab=dynamic_cast<const PGData*>(this);
    assert(ab);
    int N=ab->size1();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=ab->radials1[ia]->Integrate(t2C,ab->radials1[ib],ab->pols1[ia],ab->pols1[ib],cache,cl)*ab->ns1[ia]*ab->ns1[ib];

    return s;
}

ERI3<double> Orbital_IE::MakeOverlap3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const PGData*>(&_c);
    int Nc=c->size1();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        rsmat_t s=Integrate(qchem::Overlap3C,c->radials1[ic],c->pols1[ic]);
        s*=c->ns1[ic];
        s3.push_back(s);
    } 
    return s3;   
}
ERI3<double> Orbital_IE::MakeRepulsion3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const PGData*>(&_c);
    int Nc=c->size1();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        rsmat_t s=Integrate(qchem::Repulsion3C,c->radials1[ic],c->pols1[ic]);
        s*=c->ns1[ic];
        s3.push_back(s);
    }    
    return s3;
}
rsmat_t Orbital_IE::Integrate(qchem::IType3C type , const RadialFunction* rc, const Polarization& pc) const
{
    auto ab=dynamic_cast<const PGData*>(this);
    int N=ab->size1();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=rc->Integrate(type,ab->radials1[ia],ab->radials1[ib],ab->pols1[ia],ab->pols1[ib],pc,cache)*ab->ns1[ia]*ab->ns1[ib];
        
    return s;    
}



ERI4 Orbital_IE::MakeDirect  (const obs_t& _c) const
{
    const PGData* a=dynamic_cast<const PGData* >(this);
    const PGData* c=dynamic_cast<const PGData* >(&_c);
    assert(a);
    assert(c);
    size_t Na=a->size1(), Nc=c->size1();
    ERI4 J(Na,Nc);
    
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(ia,Na))
        {
            rsmat_t& Jab=J(ia,ib);
            for (size_t ic:iv_t(0,Nc))
                for (size_t id:iv_t(ic,Nc))
                {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=a->ns1[ia]*a->ns1[ib]*c->ns1[ic]*c->ns1[id];
                        assert(c->radials1[id]);
                        Jab(ic,id)=norm * c->radials1[id]->Integrate(a->radials1[ia],a->radials1[ib],c->radials1[ic],a->pols1[ia],a->pols1[ib],c->pols1[ic],c->pols1[id],cache);
                }
        }
    return J;
}

ERI4 Orbital_IE::MakeExchange(const obs_t& _b) const
{
    const PGData* a=dynamic_cast<const PGData* >(this);
    const PGData* b=dynamic_cast<const PGData* >(&_b);
    assert(a);
    assert(b);
    size_t Na=a->size1(), Nb=b->size1();
    ERI4 K(Na,Nb);
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(0,Nb))
           
            for (size_t ic:iv_t(ia,Na))
            {
                rsmat_t& Kac=K(ia,ic);
                for (size_t id:iv_t(0,Nb))
                {
                  //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                    double norm=a->ns1[ia]*b->ns1[ib]*a->ns1[ic]*b->ns1[id];
                    assert(b->radials1[id]);
                    if (ib==id)
                        Kac(ib,id)=norm * b->radials1[id]->Integrate(a->radials1[ia],b->radials1[ib],a->radials1[ic],a->pols1[ia],b->pols1[ib],a->pols1[ic],b->pols1[id],cache);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials1[id]->Integrate(a->radials1[ia],b->radials1[ib],a->radials1[ic],a->pols1[ia],b->pols1[ib],a->pols1[ic],b->pols1[id],cache);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials1[id]->Integrate(a->radials1[ia],b->radials1[ib],a->radials1[ic],a->pols1[ia],b->pols1[ib],a->pols1[ic],b->pols1[id],cache);
                }        
            }
    return K;
}

} //namespace PolarizedGaussian
