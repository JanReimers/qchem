// File: AtomIEClient.C Common IE client code for all atom basis sets and IEs.

#include "Imp/BasisSet/Atom/IEC.H"


void AtomIrrepIEClient::Init(const Vector<double>& exponents,const Vector<double>& norms,size_t _l)
{
    n=_l+1;
    l=_l;

    es=exponents;
    ns=norms;
}

void AtomIrrepIEClient::Init(const Vector<double>& exponents,const Vector<double>& norms,size_t _l, const std::vector<int>& _ml)
{
    n=_l+1;
    l=_l;
    ml=_ml;
    assert(ml.size()>0);

    es=exponents;
    ns=norms;
}

AtomIrrepIEClient* AtomIrrepIEClient::dcast(::IrrepIEClient* iec)
{
    assert(iec);
    AtomIrrepIEClient* aiec=dynamic_cast< AtomIrrepIEClient*>(iec);
    assert(aiec);
    return aiec;
}
const AtomIrrepIEClient* AtomIrrepIEClient::dcast(const ::IrrepIEClient* iec)
{
    assert(iec);
    const AtomIrrepIEClient* aiec=dynamic_cast< const AtomIrrepIEClient*>(iec);
    assert(aiec);
    return aiec;
}
