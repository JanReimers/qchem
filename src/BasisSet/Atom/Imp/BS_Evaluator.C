// File: BasisSet/Atom/Imp/BS_Evaluator.C  Generic hot loop 2 electron repulsion integrals (2ERIs).
module;
#include <vector>
#include <cassert>
#include <iostream>
#include <blaze/Math.h>
module qchem.BasisSet.Atom.Gaussian.NR.BS_Evaluator;
import qchem.BasisSet.Atom.Internal.AngularIntegrals; 

BS_Evaluator::rvec11_t BS_Evaluator::Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const
{
    rvec11_t Ak(0.0);
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
        Ak/=(double)nac;
    }
    return Ak;
    
}

BS_Evaluator::rvec11_t BS_Evaluator::ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* b) const
{
    rvec11_t Ak(0.0);
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
    assert(a);
    assert(c);
    size_t spanab=a->maxSpan(),spancd=c->maxSpan();
    size_t Na=a->size(), Nc=c->size();
    int la=a->Getl(), lc=c->Getl();
    rvec11_t Akac=Coulomb_AngularIntegrals(a,c);
    ERI4 J(Na,Nc);
    rvec_t na=a->Norm(), nc=c->Norm();

    for (size_t ia:a->indices())
    {
        loop_1(a->es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c->indices())
        {
            loop_2(c->es_index(ic));
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                if (ib>ia+spanab) continue;
                rsmat_t& Jab=J(ia,ib);
                loop_3(a->es_index(ib));
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    if (id>ic+spancd) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        std::cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        std::cout << Jab(ic,id) << std::endl;    
                        assert(false);
                    }
                    double norm=na[ia]*na[ib]*nc[ic]*nc[id];
                    double AkacRkac=loop_4_direct(c->es_index(id),la,lc,Akac);
                    Jab(ic,id)=AkacRkac*norm;
                }
            }
        }
    }
    return J;
}
ERI4 BS_Evaluator::Exchange(const IBS_Evaluator* a, const IBS_Evaluator* c) const
{
    assert(a);
    assert(c);
    size_t spanab=a->maxSpan(),spancd=c->maxSpan();
    size_t Na=a->size(), Nc=c->size();
    ERI4 K(Na,Nc);
    rvec_t na=a->Norm(), nc=c->Norm();
    int la=a->Getl(), lc=c->Getl();
    rvec11_t Akac=ExchangeAngularIntegrals(a,c);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c->indices())
        {
            if (ia>ic+spanab || ic>ia+spancd) continue;
            for (size_t ib:a->indices(ia))
            {
                rsmat_t& Kab=K(ia,ib);
                loop_2(a->es_index(ib));
                loop_3(c->es_index(ic));
                for (size_t id:c->indices())
                {
                    if (id>ib+spancd || ib>id+spanab) continue;
                    double norm=na[ia]*na[ib]*nc[ic]*nc[id]; 
                    double AkacRKac=loop_4_exchange(c->es_index(id),la,lc,Akac);
                    if (ic==id)
                        Kab(ic,id)=AkacRKac*norm; 
                    else if (id<ic)
                        Kab(id,ic)+=0.5*AkacRKac*norm; 
                    else
                        Kab(ic,id)+=0.5*AkacRKac*norm; 

                }
            }
        }
    }

    return K;

}
