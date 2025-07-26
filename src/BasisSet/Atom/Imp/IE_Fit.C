// File: AtomIE.C Common Fit IE code for all atom basis sets.
module;
#include <cassert>
#include <memory>
module qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;

Vector<double>  AtomIE_Fit::MakeCharge() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this); //Cross cast
    assert(a);
    Vector<double>  c(a->size());
    for (auto i:a->es.indices())  c(i)=Charge(a->es(i),a->l)*a->ns(i);
    return c;
}
SMatrix<double> AtomIE_Fit::MakeRepulsion() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this); //Cross cast
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Repulsion(a->es(i),a->es(j),l,l)*a->ns(i)*a->ns(j);

    return H;
}
Matrix<double>   AtomIE_Fit::MakeRepulsion(const Fit_IBS& _b) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);  //Cross cast
    const AtomIrrepIEClient* b=AtomIrrepIEClient::dcast(&_b);
    assert(a);
    size_t Na=a->es.size(), Nb=b->es.size();
    Matrix<double> s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=this->Repulsion(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*a->ns(j);

    return s;
}

