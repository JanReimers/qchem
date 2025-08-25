// File: BasisSet/Atom/BS_Evaluator.C Create Rk structures for 2 electron repulsion integrals (2ERIs).
module;
#include <vector>
export module BasisSet.Atom.BS_Evaluator;
export import BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Atom.Internal.AngularIntegrals;


export class BS_Evaluator : public Cache4
{
public:
    virtual ~BS_Evaluator() {};
    virtual void Register(IBS_Evaluator*)=0;
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