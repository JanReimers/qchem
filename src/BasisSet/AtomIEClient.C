// File: AtomIEClient.C Common IE client code for all atom basis sets and IEs.

#include "Imp/BasisSet/AtomIEClient.H"

template <class T> void FillPower(Vector<T>& arr,T start, T stop);

void AtomIrrepIEClient::Init(double minexp,double maxexp,size_t _l, int _m)
{
    n=_l+1;
    l=_l;
    m=_m;
    
    FillPower(es,minexp,maxexp);
    for (auto i:es.indices())  ns(i)=Norm(es(i),l);
}

void AtomIrrepIEClient::Init(std::set<double>& exponents,size_t _l, int _m)
{
    n=_l+1;
    l=_l;
    m=_m;

    int i=1;
    for (auto& e:exponents) es(i++)=e;
    for (auto i:es.indices())  ns(i)=Norm(es(i),l);
}

size_t AtomIEClient::LMax() const 
{
    return Max(Ls);
}
void AtomIEClient::Append(const AtomIrrepIEClient* ic)
{
    itsIrreps.push_back(ic);
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ls.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ls(j)=ic->l;
        BFGrouper::Append(ic->es(i),ic->l,j);
    }

}

