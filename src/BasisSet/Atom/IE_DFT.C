// File: Atom/IE_DFT.C Common DFT IE code for all atom basis sets.

#include "Imp/BasisSet/Atom/IE_DFT.H"
#include "Imp/BasisSet/Atom/IEC.H"
#include <Fit_IBS.H>

template <class T> typename AtomIE_DFT<T>::ERI3 AtomIE_DFT<T>::MakeOverlap3C  (const fbs_t& _c) const
{
    auto& c=AtomIrrepIEClient::dcast(_c);
    ERI3 s3;
    for (auto i:c.indices()) s3.push_back(MakeOverlap(c(i)));
    return s3;
}
template <class T> typename AtomIE_DFT<T>::ERI3 AtomIE_DFT<T>::MakeRepulsion3C(const fbs_t& _c) const
{
    auto& c=AtomIrrepIEClient::dcast(_c);
    ERI3 s3;
    for (auto i:c.indices()) s3.push_back(MakeRepulsion(c(i)));
    return s3;
}
template <class T> typename AtomIE_DFT<T>::SMat AtomIE_DFT<T>::MakeOverlap  (const bf_tuple& c) const
{    
    const AtomIrrepIEClient* ab=dynamic_cast<const AtomIrrepIEClient*>(this); //cross cast.
    assert(ab);
    size_t N=ab->size();
    int Nc,Lc;
    double ec,nc;
    std::tie(Nc,Lc,ec,nc)=c;
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=this->Overlap(ab->es(i)+ab->es(j),ec,ab->l+ab->l+Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}
template <class T> typename AtomIE_DFT<T>::SMat AtomIE_DFT<T>::MakeRepulsion(const bf_tuple& c) const
{    
    const AtomIrrepIEClient* ab=dynamic_cast<const AtomIrrepIEClient*>(this); //cross cast.
    assert(ab);
    size_t N=ab->size();
    int Nc,Lc;
    double ec,nc;
    std::tie(Nc,Lc,ec,nc)=c;
    SMat s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=this->Repulsion(ab->es(i)+ab->es(j),ec,ab->l,Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}

template class AtomIE_DFT<double>;

