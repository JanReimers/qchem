// File: AtomIEClient.C Common IE client code for all atom basis sets and IEs.

#include "Imp/BasisSet/AtomIEClient.H"

template <class T> void FillPower(Vector<T>& arr,T start, T stop);

void AtomIrrepIEClient::Init(double minexp,double maxexp,size_t L, int m)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      Fill(Ms,m);
      for (auto i:es.indices())  ns(i)=Norm(es(i),Ls(i));
}

void AtomIrrepIEClient::Init(std::set<double>& exponents,size_t L, int m)
{
    int i=1;
    for (auto& e:exponents) es(i++)=e;
    Fill(Ns,L+1);
    Fill(Ls,L);
    Fill(Ms,m);
    for (auto i:es.indices())  ns(i)=Norm(es(i),Ls(i));
}

void AtomIEClient::Append(const AtomIrrepIEClient* ic)
{
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ns.SetLimits(N,true);
    Ls.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ns(j)=ic->Ns(i);
        Ls(j)=ic->Ls(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
        BFGrouper::Append(es(j),Ls(j),j);
    }

}

