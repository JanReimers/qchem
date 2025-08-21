// File: AtomIE.C Common IE code for all atom basis sets.
module;
#include <vector>
#include <iostream>
#include <cassert>
module qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;

template <class T> SMatrix<T> AtomIE_Nuclear <T>::MakeNuclear(const Cluster* cl) const
{
    assert(cl);
    assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
    int Z=-cl->GetNuclearCharge(); 
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);  // Cross cast
    assert(a);
 
    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Z*pie->Inv_r1(a->es(i),a->es(j),2*l)*a->ns(i)*a->ns(j);

    return H;
}

template class AtomIE_Nuclear<double>;

