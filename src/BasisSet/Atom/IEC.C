// File: AtomIEClient.C Common IE client code for all atom basis sets and IEs.

#include "Imp/BasisSet/Atom/IEC.H"


void AtomIrrepIEClient::Init(const Vector<double>& exponents,const Vector<double>& norms,size_t _l, int _m)
{
    n=_l+1;
    l=_l;
    m=_m;

    es=exponents;
    ns=norms;
}

