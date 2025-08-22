// File: AtomIE.C Common Fit IE code for all atom basis sets.
module;
#include <cassert>
#include <memory>
module qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;

Matrix<double>   AtomIE_Fit::MakeRepulsion(const Fit_IBS& _b) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);  //Cross cast
    const AtomIrrepIEClient* b=AtomIrrepIEClient::dcast(&_b);
    assert(a);
    size_t Na=a->es.size(), Nb=b->es.size();
    Matrix<double> s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=pie->Repulsion(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*a->ns(j);

    return s;
}

