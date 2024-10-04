#include "ERIProxy.H"

template <class T> TERIProxy<T>::TERIProxy(T& eril,int start_ab, int start_cd)
 : itsMasterERIList(eril)
 , itsStart_a(start_ab)
 , itsStart_b(start_ab)
 , itsStart_c(start_cd)
 , itsStart_d(start_cd)
{
    assert(itsStart_a>0);
    assert(itsStart_b>0);
    assert(itsStart_c>0);
    assert(itsStart_d>0);
    //ctor
}

template <class T> TERIProxy<T>::TERIProxy(T& eril,int sa,int sb,int sc, int sd)
 : itsMasterERIList(eril)
 , itsStart_a(sa)
 , itsStart_b(sb)
 , itsStart_c(sc)
 , itsStart_d(sd)
{
    assert(itsStart_a>0);
    assert(itsStart_b>0);
    assert(itsStart_c>0);
    assert(itsStart_d>0);
    //ctor
}

template class TERIProxy<ERIList>;
template class TERIProxy<ERIList1>;
