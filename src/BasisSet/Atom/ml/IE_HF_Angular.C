// File: Atom_ml_IE_HF_Angular.H  Angular 2e-Integrals for atom-ml HF basis sets.

#include "Imp/BasisSet/Atom/ml/IE_HF_Angular.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/BasisSet/Atom/IEC.H"

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

IE_BS_2E_Angular::RVec IE_BS_2E_Angular::Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const
 {
    if (a->ml.size()==0 && c->ml.size()==0)
        return AngularIntegrals::Coulomb(a->l,c->l,a->m,c->m); //old style
    else
    {
        RVec Ak;
        size_t nac=a->ml.size()*c->ml.size();
        
        for (auto ma:a->ml)
        for (auto mc:c->ml)
            Ak+=AngularIntegrals::Coulomb(a->l,c->l,ma,mc);
        return Ak/nac;
    }
 }
IE_BS_2E_Angular::RVec IE_BS_2E_Angular::ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const
{
    if (a->ml.size()==0 && b->ml.size()==0)
        return AngularIntegrals::Exchange(a->l,b->l,a->m,b->m);
    else
    {
        RVec Ak;
        size_t nab=(a->ml.size())*(b->ml.size());
        for (auto ma:a->ml)
        for (auto mb:b->ml)
            Ak+=AngularIntegrals::Exchange(a->l,b->l,ma,mb);
        return Ak/nab;
    }
}

}

