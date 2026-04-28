module;
#include <cassert>
#include <vector>
module qchem.BasisSet.Internal.DB_Cache;



template <class T> void DB_BS_HF<T>::Append(const Orbital_HF_IBS<T>* oibs)
{
    itsIrreps.push_back(oibs);
}
template <class T> ERI4 DB_BS_HF<T>::Direct(IDType a,IDType c) const
{
    // assert(a<=c);
    if (Jac.size()==0) MakeDirect();
    //cout << "GetRepulsion4C_new a,c=" << a.GetIndex() << " " << c.GetIndex() << endl;
    assert(Jac.find(a)!=Jac.end());
    assert(Jac[a].find(c)!=Jac[a].end());
    
    return Jac[a][c];
}
template <class T> ERI4 DB_BS_HF<T>::Exchange(IDType a,IDType b) const
{
    // assert(a<=b);
    if (Kab.size()==0) MakeExchange(); 
    //cout << "GetExchange4C_new a,b=" << a.GetIndex() << " " << b.GetIndex() << endl;
    assert(Kab.find(a)!=Kab.end());
    assert(Kab[a].find(b)!=Kab[a].end());
    
    return Kab[a][b];
}
template <class T> void DB_BS_HF<T>::MakeDirect() const
{
    Jac.clear();
    for (auto a: itsIrreps)
        for (auto c: itsIrreps) //TODO run from ia n
        {
            IDType ia=a->GetID(), ic=c->GetID();
            ERI4 jac=MakeDirect(a,c);
            Jac[ia][ic]=jac;
            Jac[ic][ia]= ia==ic ? jac : jac.Transpose();
        }
}
template <class T> void DB_BS_HF<T>::MakeExchange() const
{
    Kab.clear();
    for (auto a: itsIrreps)
        for (auto b: itsIrreps) 
        {
            IDType ia=a->GetID(), ib=b->GetID();
            ERI4 kab=MakeExchange(a,b);
            Kab[ia][ib]=kab;
            Kab[ib][ia]= ia==ib ? kab : kab.Transpose();
        }
    
}

template class DB_BS_HF<double>;
