// File: AtomIEClient.C Common IE client code for all atom basis sets and IEs.

#include "Imp/BasisSet/AtomIEClient.H"

template <class T> void FillPower(Vector<T>& arr,T start, T stop);

void AtomIEClient::Init(double minexp,double maxexp,size_t L, int m)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      Fill(Ms,m);
      for (auto i:es.indices())  ns(i)=Norm(es(i),Ls(i));
}

void AtomIEClient::Init(std::set<double>& exponents,size_t L, int m)
{
    int i=1;
    for (auto& e:exponents) es(i++)=e;
    Fill(Ns,L+1);
    Fill(Ls,L);
    Fill(Ms,m);
    for (auto i:es.indices())  ns(i)=Norm(es(i),Ls(i));
}


