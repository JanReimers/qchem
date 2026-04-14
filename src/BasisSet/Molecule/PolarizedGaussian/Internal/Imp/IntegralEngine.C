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
    for (auto r:a->radials)
    {
        c[i]=r->GetCharge(a->pols[i])*a->ns[i]; 
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
            s(ia,ib)=a->radials[ia]->Integrate(Repulsion2C,
                b->radials[ib],a->pols[ia],b->pols[ib],cache)*a->ns[ia]*b->ns[ib];
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
            s(ia,ib)=ab->radials[ia]->Integrate(t2C,ab->radials[ib],ab->pols[ia],ab->pols[ib],cache,cl)*ab->ns[ia]*ab->ns[ib];

    return s;
}

ERI3<double> Orbital_IE::MakeOverlap3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const PGData*>(&_c);
    int Nc=c->size1();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        rsmat_t s=Integrate(qchem::Overlap3C,c->radials[ic],c->pols[ic]);
        s*=c->ns[ic];
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
        rsmat_t s=Integrate(qchem::Repulsion3C,c->radials[ic],c->pols[ic]);
        s*=c->ns[ic];
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
            s(ia,ib)=rc->Integrate(type,ab->radials[ia],ab->radials[ib],ab->pols[ia],ab->pols[ib],pc,cache)*ab->ns[ia]*ab->ns[ib];
        
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
                        double norm=a->ns[ia]*a->ns[ib]*c->ns[ic]*c->ns[id];
                        assert(c->radials[id]);
                        Jab(ic,id)=norm * c->radials[id]->Integrate(a->radials[ia],a->radials[ib],c->radials[ic],a->pols[ia],a->pols[ib],c->pols[ic],c->pols[id],cache);
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
                    double norm=a->ns[ia]*b->ns[ib]*a->ns[ic]*b->ns[id];
                    assert(b->radials[id]);
                    if (ib==id)
                        Kac(ib,id)=norm * b->radials[id]->Integrate(a->radials[ia],b->radials[ib],a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id],cache);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials[id]->Integrate(a->radials[ia],b->radials[ib],a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id],cache);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials[id]->Integrate(a->radials[ia],b->radials[ib],a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id],cache);
                }        
            }
    return K;
}

} //namespace PolarizedGaussian
