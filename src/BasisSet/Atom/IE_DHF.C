// File: AtomIE.C Common DHF IE code for all atom basis sets.

#include "Imp/BasisSet/Atom/IE_DHF.H"
#include "Imp/BasisSet/Atom/IEC.H"
#include <DHF_IBS.H>

template <class T> typename Integrals_Base<T>::Mat  AtomIE_XKinetic<T>::MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    const AtomIrrepIEClient* b=dynamic_cast<const AtomIrrepIEClient*>(rkbs);
    assert(a->l==b->l);
    size_t l=a->l;
    size_t Na=a->size();
    size_t Nb=b->size();
    Matrix<double> Hk(Na,Nb);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols())
            Hk(i,j)=(Grad2(a->es(i),b->es(j),l,l) + l*(l+1)*Inv_r2(a->es(i),b->es(j),2*l))*a->ns(i)*b->ns(j);

    return Hk;
}

template class AtomIE_RKBL<double>;
template class AtomIE_RKBS<double>;
