// File: BasisSet/Atom/BS_Evaluator.C Create Rk structures for 2 electron repulsion integrals (2ERIs).
module;
#include <vector>
#include <cassert>
#include <iostream>
export module BasisSet.Atom.BS_Evaluator;
export import BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Atom.Internal.AngularIntegrals;
import qchem.BasisSet.Internal.ERI4;


export class BS_Evaluator : public Cache4
{
public:
    virtual ~BS_Evaluator() {};
    virtual void Register(IBS_Evaluator*)=0;
    ERI4 Direct  (const IBS_Evaluator* a, const IBS_Evaluator* c) const;
    ERI4 Exchange(const IBS_Evaluator* a, const IBS_Evaluator* c) const;
    virtual RVec Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    virtual RVec ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    virtual RVec loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual RVec loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
};

// template <class T> const Vector<T>& operator+=(Vector<T>& a, const Vector<T>& b);
// {
//     if (a.size()==0) 
//     {
//         a.SetLimits(b.GetLimits());
//         Fill(a,0.0);
//     }
//     return ArrayAdd(a,b);
// }

RVec BS_Evaluator::Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const
{
    RVec Ak;
    int la=a->Getl(),lc=c->Getl();
    const IBS_Evaluator::is_t& amls=a->Getmls(),cmls=c->Getmls();
    size_t nac=amls.size()*cmls.size();
    size_t g=(2*la+1)*(2*lc+1); //degenracy
    if (nac==g || nac==0)
    {
        Ak+=AngularIntegrals::Coulomb(la,lc);
    }
    else
    {
        // For direct integrals these actually factor.  But for exchange they do not.
        // So it may not be worth while coding the factored version here.
        for (auto ma:amls)
        for (auto mc:cmls)
            Ak+=AngularIntegrals::Coulomb(la,lc,ma,mc);
        Ak/=nac;
    }
    return Ak;
    
}

RVec BS_Evaluator::ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* b) const
{
    RVec Ak;
    int la=a->Getl(),lb=b->Getl();
    const IBS_Evaluator::is_t& amls=a->Getmls(),bmls=b->Getmls();
    size_t nab=amls.size()*bmls.size();
    size_t g=(2*la+1)*(2*lb+1); //degenracy
    if (nab==0 || nab==g)
    {
        Ak+=AngularIntegrals::Exchange(la,lb);
    }
    else
    {
        for (auto ma:amls)
        for (auto mb:bmls)
            Ak+=AngularIntegrals::Exchange(la,lb,ma,mb);
        Ak/=nab;
    }
    return Ak;
}

ERI4 BS_Evaluator::Direct  (const IBS_Evaluator* a, const IBS_Evaluator* c) const
{
    using SMat=IBS_Evaluator::SMat;
    using ds_t=IBS_Evaluator::ds_t;
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    ds_t na=a->Norm(), nc=c->Norm();

    for (size_t ia:a->indices())
    {
        loop_1(a->es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c->indices())
        {
            loop_2(c->es_index(ic));
            int la=a->Getl(), lc=c->Getl();
            RVec Akac=Coulomb_AngularIntegrals(a,c);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia+1,ib+1);
                loop_3(a->es_index(ib));
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    if (Jab(ic+1,id+1)!=0.0)
                    {
                        std::cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        std::cout << Jab(ic+1,id+1) << std::endl;    
                        assert(false);
                    }
                    double norm=na[ia]*na[ib]*nc[ic]*nc[id];
                    RVec Rkac=loop_4_direct(c->es_index(id),la,lc);
                    Jab(ic+1,id+1)=Akac*Rkac*norm;
                }
            }
        }
    }
    return J;
}
ERI4 BS_Evaluator::Exchange(const IBS_Evaluator* a, const IBS_Evaluator* c) const
{
    using SMat=IBS_Evaluator::SMat;
    using ds_t=IBS_Evaluator::ds_t;
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 K(Na,Nc);
    ds_t na=a->Norm(), nc=c->Norm();
    for (size_t ia:a->indices())
    {
        loop_1(a->es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c->indices())
        {
            int la=a->Getl(), lc=c->Getl();
            RVec Akac=ExchangeAngularIntegrals(a,c);
            for (size_t ib:a->indices(ia))
            {
                SMat& Kab=K(ia+1,ib+1);
                loop_2(a->es_index(ib));
                loop_3(c->es_index(ic));
                for (size_t id:c->indices())
                {
                    double norm=na[ia]*na[ib]*nc[ic]*nc[id]; 
                    RVec RKac=loop_4_exchange(c->es_index(id),la,lc);
                    if (ic==id)
                        Kab(ic+1,id+1)=Akac*RKac*norm; 
                    else if (id<ic)
                        Kab(id+1,ic+1)+=0.5*Akac*RKac*norm; 
                    else
                        Kab(ic+1,id+1)+=0.5*Akac*RKac*norm; 

                }
            }
        }
    }

    return K;

}
