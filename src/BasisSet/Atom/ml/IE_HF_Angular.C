// File: Atom_ml_IE_HF_Angular.C  Angular 2e-Integrals for atom-ml HF basis sets.
module;
#include <vector>
export module qchem.BasisSet.Atom.Internal.ml.Angular;

import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.Internal.AngularIntegrals;
export import qchem.Types;

template <class T> const Vector<T>& operator+=(Vector<T>& a, const Vector<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    return ArrayAdd(a,b);
}

namespace Atom_ml
{
export class IE_BS_2E_Angular : public virtual ::AtomIE_BS_2E_Angular
{
public:
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const;
};

RVec IE_BS_2E_Angular::Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const
 {
    RVec Ak;
    size_t nac=a->ml.size()*c->ml.size();
    size_t g=(2*a->l+1)*(2*c->l+1); //degenracy
    if (nac==g)
    {
        Ak+=AngularIntegrals::Coulomb(a->l,c->l);
    }
    else
    {
        // For direct integrals these actually factor.  But for exchange they do not.
        // So it may not be worth while coding the factored version here.
        for (auto ma:a->ml)
        for (auto mc:c->ml)
            Ak+=AngularIntegrals::Coulomb(a->l,c->l,ma,mc);
        Ak/=nac;
    }
    return Ak;
    
 }
RVec IE_BS_2E_Angular::ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const
{
    
    RVec Ak;
    size_t nab=(a->ml.size())*(b->ml.size());
    size_t g=(2*a->l+1)*(2*b->l+1); //degenracy
    if (nab==g)
    {
        Ak+=AngularIntegrals::Exchange(a->l,b->l);
    }
    else
    {
        for (auto ma:a->ml)
        for (auto mb:b->ml)
            Ak+=AngularIntegrals::Exchange(a->l,b->l,ma,mb);
        Ak/=nab;
    }

    return Ak;
    
}

}

