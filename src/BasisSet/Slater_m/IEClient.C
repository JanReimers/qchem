
#include "Imp/BasisSet/Slater_m/IEClient.H"
#include "Imp/Integrals/SlaterIntegrals.H"

using std::cout;
using std::endl;

template <class T> void FillPower(Vector<T>& arr,T start, T stop);
//template void FillPower(Vector<double>& arr,double start, double stop);

namespace Slater_m
{
    
void IrrepIEClient::Init(double minexp,double maxexp,size_t L, int m)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      Fill(Ms,m);
      for (auto i:es.indices())  ns(i)=SlaterNorm(es(i),Ns(i));
}

void IEClient::Append(const IrrepIEClient* ic)
{
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ns.SetLimits(N,true);
    Ls.SetLimits(N,true);
    Ms.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ns(j)=ic->Ns(i);
        Ls(j)=ic->Ls(i);
        Ms(j)=ic->Ms(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
        
        size_t index=unique_es.size();
        if (const auto &ie =unique_es.find(es(j));ie==unique_es.end())
        {
            unique_esv.push_back(es(j));
            unique_es[es(j)]=index;
        }
        else 
            index=ie->second;
        
        es_indices.push_back(index);
        
        if (const auto &il =L_indices.find(Ls(j));il==L_indices.end())
            L_indices[Ls(j)]=std::vector<size_t>();
        
        L_indices[Ls(j)].push_back(j);
        
        
    }
//    for (auto e:unique_esv) cout << e << " ";
//    cout << endl;
//    for (auto e:unique_es) cout << e.first << " ";
//    cout << endl;
//    for (auto i:es_indexes) cout << i << " ";
//    cout << endl;

}

const std::vector<size_t>& IEClient::indices(size_t l) const
{
    auto i=L_indices.find(l);
    assert(i!=L_indices.end());
    return i->second;
}

using std::cout;
using std::endl;


const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}




} //namespace
